"""
pn1d_reference.py -- reference numbers + plots for the 1D p+/n junction.

Same parameters as the COMSOL model (comsol/pn1d.mph). Does three things:
  1. depletion-approximation numbers (Vbi, W, xn, xp, Emax, J0 ...)
  2. an independent Newton solve of Poisson at equilibrium (Boltzmann stats),
     finite volumes on the same kind of graded mesh -> something to diff COMSOL against
  3. the I-V I expect (Shockley short diode + SRH recombination/generation)

usage:  python analytics/pn1d_reference.py
writes: docs/img/*.png, data/reference/*.json|csv
"""

from pathlib import Path
import json

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
FIG = ROOT / "docs" / "img"
DATA = ROOT / "data" / "reference"
FIG.mkdir(parents=True, exist_ok=True)
DATA.mkdir(parents=True, exist_ok=True)

# --- constants (SI)
q = 1.602176634e-19
kB = 1.380649e-23
eps0 = 8.8541878128e-12

# --- model parameters (keep in sync with COMSOL Parameters node!)
T0 = 300.0
Vt = kB * T0 / q
NA = 1e19 * 1e6          # m^-3
ND = 1e16 * 1e6
ni = 9.65e9 * 1e6
Lp_reg = 2e-6
Ln_reg = 10e-6
L_tot = Lp_reg + Ln_reg
A_d = 1e-4 * 1e-4        # m^2
eps_s = 11.7 * eps0
mun = 1400e-4            # m^2/(V s), constant like in COMSOL (for now)
mup = 470e-4
tau_n = 1e-6
tau_p = 1e-6

# minority mobilities for J0 only (doping dependent, Masetti). COMSOL still uses constants
mup_n_side = 430e-4      # lacune nel lato n (1e16 cm^-3)
mun_p_side = 115e-4      # elettroni nel lato p+ (1e19 cm^-3)


def analytic(V=0.0):
    """Depletion approx at anode bias V."""
    Vbi = Vt * np.log(NA * ND / ni**2)
    W = np.sqrt(2 * eps_s * (Vbi - V) / q * (1 / NA + 1 / ND))
    xn = W * NA / (NA + ND)
    xp = W * ND / (NA + ND)
    Emax = 2 * (Vbi - V) / W
    return dict(Vbi=Vbi, W=W, xn=xn, xp=xp, Emax=Emax)


def saturation_currents():
    Dp = mup_n_side * Vt
    Dn = mun_p_side * Vt
    Lp = np.sqrt(Dp * tau_p)
    Ln = np.sqrt(Dn * tau_n)
    a0 = analytic(0.0)
    Wn = Ln_reg - a0["xn"]          # neutral n region
    Wp = Lp_reg - a0["xp"]          # neutral p+ region
    Jp0 = q * ni**2 / ND * Dp / Lp / np.tanh(Wn / Lp)
    Jn0 = q * ni**2 / NA * Dn / Ln / np.tanh(Wp / Ln)
    return dict(Dp=Dp, Dn=Dn, Lp=Lp, Ln=Ln, Wn=Wn, Jp0=Jp0, Jn0=Jn0, J0=Jp0 + Jn0)


# --- numerical Poisson
def solve_poisson_equilibrium(n_left=1500, n_right=2500):
    """-d/dx(eps dpsi/dx) = q(p - n + Ndop), phin = phip = 0."""
    # geometric grid, fine at the junction on both sides (same idea as the COMSOL Distribution nodes)
    def geo(a, b, n, ratio, reverse):
        r = ratio ** (1 / (n - 1))
        h = np.array([r**k for k in range(n)])
        h = h / h.sum() * (b - a)
        if reverse:
            h = h[::-1]
        return a + np.concatenate(([0.0], np.cumsum(h)))

    xl = geo(0.0, Lp_reg, n_left, 200, reverse=True)
    xr = geo(Lp_reg, L_tot, n_right, 1000, reverse=False)
    xl[-1] = Lp_reg
    x = np.concatenate((xl, xr[1:]))
    i_junction = len(xl) - 1
    Ndop = np.where(np.arange(len(x)) < i_junction, -NA, ND).astype(float)
    Ndop[i_junction] = 0.5 * (ND - NA)

    psi = Vt * np.arcsinh(Ndop / (2 * ni))
    psi_L = Vt * np.arcsinh(-NA / (2 * ni))
    psi_R = Vt * np.arcsinh(ND / (2 * ni))
    psi[0], psi[-1] = psi_L, psi_R

    h = np.diff(x)
    hm = 0.5 * (h[:-1] + h[1:])           # control volume of interior nodes
    for it in range(200):
        n = ni * np.exp(psi / Vt)
        p = ni * np.exp(-psi / Vt)
        flux = eps_s * np.diff(psi) / h
        F = (flux[1:] - flux[:-1]) + q * (p - n + Ndop)[1:-1] * hm
        dF_dpsi = -q * (p + n)[1:-1] / Vt * hm
        main = -eps_s / h[1:] - eps_s / h[:-1] + dF_dpsi
        lower = eps_s / h[1:-1]
        upper = eps_s / h[1:-1]
        # tridiagonal solve (Thomas). could use scipy but whatever
        nN = len(main)
        c = np.zeros(nN)
        d = np.zeros(nN)
        rhs = -F
        c[0] = upper[0] / main[0]
        d[0] = rhs[0] / main[0]
        for i in range(1, nN):
            den = main[i] - lower[i - 1] * c[i - 1]
            c[i] = upper[i] / den if i < nN - 1 else 0.0
            d[i] = (rhs[i] - lower[i - 1] * d[i - 1]) / den
        dpsi = np.zeros(nN)
        dpsi[-1] = d[-1]
        for i in range(nN - 2, -1, -1):
            dpsi[i] = d[i] - c[i] * dpsi[i + 1]
        dpsi = np.clip(dpsi, -5 * Vt, 5 * Vt)  # damped Newton, otherwise first iterations blow up
        psi[1:-1] += dpsi
        if np.max(np.abs(dpsi)) < 1e-10:
            break
    n = ni * np.exp(psi / Vt)
    p = ni * np.exp(-psi / Vt)
    xm = 0.5 * (x[1:] + x[:-1])
    E = -np.diff(psi) / h
    return dict(x=x, psi=psi, n=n, p=p, Ndop=Ndop, xm=xm, E=E, iterations=it + 1)


# --- expected I-V
def iv_curves(J0):
    V = np.linspace(0.0, 0.8, 161)
    out = dict(V=V, Jdiff=J0 * (np.exp(V / Vt) - 1))
    Wv = np.array([analytic(v)["W"] if v < analytic(0)["Vbi"] - 0.05 else np.nan for v in V])
    out["Jrec"] = q * Wv * ni / (2 * tau_n) * (np.exp(V / (2 * Vt)) - 1)
    Vr = np.linspace(0, -40, 161)
    Wr = np.array([analytic(v)["W"] for v in Vr])
    out["Vr"] = Vr
    out["Jgen"] = q * ni * Wr / (2 * tau_n)
    out["Jdiff_rev"] = J0 * np.ones_like(Vr)
    return out


# --- plot style
plt.rcParams.update({
    "figure.dpi": 150,
    "savefig.dpi": 150,
    "font.size": 10,
    "font.family": "monospace",
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linestyle": ":",
    "axes.titlesize": 10,
})
C1, C2, C3 = "C0", "C3", "0.4"


def save(fig, name):
    fig.tight_layout()
    fig.savefig(FIG / name)
    plt.close(fig)


def main():
    a0 = analytic(0.0)
    sat = saturation_currents()
    sol = solve_poisson_equilibrium()
    x_um = sol["x"] * 1e6
    xm_um = sol["xm"] * 1e6
    zoom = (1.9, 2.6)

    # --- numeric vs analytic
    rho = q * (sol["p"] - sol["n"] + sol["Ndop"])
    i_j = np.argmin(np.abs(sol["x"] - Lp_reg))
    dpsi_num = sol["psi"][-1] - sol["psi"][0]
    Emax_num = np.max(np.abs(sol["E"]))
    # "numerical" depletion edge: where n drops below ND/2 on the n side
    mask_n = (sol["x"] > Lp_reg) & (sol["n"] < 0.5 * ND)
    W_num_n = sol["x"][mask_n].max() - Lp_reg if mask_n.any() else np.nan

    # --- figura 1: drogaggio netto
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.semilogy(x_um, np.where(sol["Ndop"] < 0, NA, np.nan) * 1e-6, color=C2, label="$N_A$ (p+ side)")
    ax.semilogy(x_um, np.where(sol["Ndop"] > 0, ND, np.nan) * 1e-6, color=C1, label="$N_D$ (n side)")
    ax.axvline(Lp_reg * 1e6, color=C3, ls="--", lw=1)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("doping (cm$^{-3}$)")
    ax.set_title("doping profile")
    ax.set_ylim(1e15, 1e20)
    ax.legend(frameon=False, loc="center right")
    save(fig, "doping.png")

    # --- figura 2: potenziale (intero e zoom)
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.2))
    for ax, xl in zip(axs, [(0, 12), zoom]):
        ax.plot(x_um, sol["psi"], color=C1)
        ax.set_xlim(*xl)
        ax.set_xlabel("x (µm)")
        ax.set_ylabel("ψ (V)")
    axs[0].set_title("psi, whole device")
    axs[1].set_title("psi, zoom on junction")
    save(fig, "psi.png")

    # --- figura 3: carica (zoom, scala limitata)
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(x_um, rho, color=C1)
    ax.set_xlim(*zoom)
    ax.set_ylim(-2000, 2000)
    ax.axhline(q * ND, color=C3, ls=":", lw=1)
    ax.text(2.42, q * ND + 120, "$+qN_D$", color=C3)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("ρ (C/m³)")
    ax.set_title("charge density (the -qN_A spike on p+ side is off-scale)")
    save(fig, "rho.png")

    # --- figura 4: campo
    fig, ax = plt.subplots(figsize=(6.4, 3.2))
    ax.plot(xm_um, np.abs(sol["E"]) * 1e-2, color=C1, label="numerical (python Newton)")
    xt = np.array([Lp_reg - a0["xp"], Lp_reg, Lp_reg + a0["xn"]]) * 1e6
    ax.plot(xt, [0, a0["Emax"] * 1e-2, 0], color=C2, ls="--", label="depletion approx")
    ax.set_xlim(*zoom)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("|E| (V/cm)")
    ax.set_title("|E| at equilibrium (E points n -> p, so E < 0)")
    ax.legend(frameon=False)
    save(fig, "efield.png")

    # --- figura 5: concentrazioni
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.semilogy(x_um, sol["n"] * 1e-6, color=C1, label="n")
    ax.semilogy(x_um, sol["p"] * 1e-6, color=C2, label="p")
    ax.set_xlim(1.5, 3.0)
    ax.set_ylim(1e0, 1e20)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("carriers (cm$^{-3}$)")
    ax.set_title("carriers at equilibrium")
    ax.legend(frameon=False)
    save(fig, "carriers.png")

    # --- figura 6: bande
    Ei = -sol["psi"]
    Eg = 1.12
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    ax.plot(x_um, Ei + Eg / 2, color=C1, label="$E_c$")
    ax.plot(x_um, Ei - Eg / 2, color=C2, label="$E_v$")
    ax.plot(x_um, Ei, color=C3, ls=":", label="$E_i$")
    ax.axhline(0, color="k", ls="--", lw=1, label="$E_F$")
    ax.set_xlim(1.5, 3.0)
    ax.set_xlabel("x (µm)")
    ax.set_ylabel("energy (eV)")
    ax.set_title("band diagram, V = 0")
    ax.legend(frameon=False, ncol=4, loc="upper right", fontsize=8)
    save(fig, "bands.png")

    # --- figura 7: I-V attesa
    iv = iv_curves(sat["J0"])
    fig, axs = plt.subplots(1, 2, figsize=(9, 3.4))
    ax = axs[0]
    ax.semilogy(iv["V"], iv["Jdiff"] * A_d, color=C1, ls="--", label="diffusion (Shockley, short diode)")
    ax.semilogy(iv["V"], iv["Jrec"] * A_d, color=C2, ls="--", label="SRH recombination in SCR")
    ax.semilogy(iv["V"], (iv["Jdiff"] + np.nan_to_num(iv["Jrec"])) * A_d, color="k", label="sum")
    ax.set_ylim(1e-14, 1e-1)
    ax.set_xlabel("V (V)")
    ax.set_ylabel("I (A)")
    ax.set_title("expected forward I-V (no HL injection, no Rs)")
    ax.legend(frameon=False, fontsize=8)
    ax = axs[1]
    ax.semilogy(-iv["Vr"], iv["Jgen"] * A_d, color=C2, label="SRH generation in SCR")
    ax.semilogy(-iv["Vr"], iv["Jdiff_rev"] * A_d, color=C1, ls="--", label="diffusion (J0*A)")
    ax.set_xlabel("|V| reverse (V)")
    ax.set_ylabel("|I| (A)")
    ax.set_title("expected reverse current")
    ax.legend(frameon=False, fontsize=8)
    save(fig, "iv_expected.png")

    # --- dump numbers
    W10 = analytic(-10)
    W40 = analytic(-40)
    values = {
        "Vt_V": Vt,
        "Vbi_V": a0["Vbi"],
        "W0_um": a0["W"] * 1e6,
        "W0_corr_um": np.sqrt(2 * eps_s * (a0["Vbi"] - 2 * Vt) / (q * ND)) * 1e6,
        "xn0_um": a0["xn"] * 1e6,
        "xp0_nm": a0["xp"] * 1e9,
        "Emax0_Vcm": a0["Emax"] * 1e-2,
        "W_m10V_um": W10["W"] * 1e6,
        "Emax_m10V_Vcm": W10["Emax"] * 1e-2,
        "W_m40V_um": W40["W"] * 1e6,
        "Cj0_pF": eps_s / a0["W"] * A_d * 1e12,
        "Cj_m10V_pF": eps_s / W10["W"] * A_d * 1e12,
        "psi_anode_V": Vt * np.arcsinh(-NA / (2 * ni)),
        "psi_cathode_V": Vt * np.arcsinh(ND / (2 * ni)),
        "Lp_um": sat["Lp"] * 1e6,
        "Jp0_Acm2": sat["Jp0"] * 1e-4,
        "Jn0_Acm2": sat["Jn0"] * 1e-4,
        "J0_Acm2": sat["J0"] * 1e-4,
        "I0_A": sat["J0"] * A_d,
        "Jgen_m10V_Acm2": q * ni * W10["W"] / (2 * tau_n) * 1e-4,
        "numeric": {
            "delta_psi_V": dpsi_num,
            "Emax_Vcm": Emax_num * 1e-2,
            "depletion_n_side_um": W_num_n * 1e6,
            "newton_iterations": sol["iterations"],
            "nodes": len(sol["x"]),
        },
    }
    (DATA / "expected_values.json").write_text(json.dumps(values, indent=2))
    np.savetxt(
        DATA / "equilibrium_python.csv",
        np.column_stack([sol["x"], sol["psi"], sol["n"], sol["p"], rho]),
        delimiter=",",
        header="x_m,psi_V,n_m-3,p_m-3,rho_Cm-3",
        comments="",
    )
    for k, v in values.items():
        print(k, v)


if __name__ == "__main__":
    main()
