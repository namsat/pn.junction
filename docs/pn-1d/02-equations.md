# 02 - equations

Drift-diffusion model: 3 PDEs, 3 unknowns. This page is about which equations go in as PDEs, which don't, and how to type them into COMSOL's *Coefficient Form PDE*.

## what goes in as a PDE and what doesn't

This confused me at first. Why not put in the equation for n, Shockley, Fermi-Dirac...? Because they're not the same kind of thing:

| kind | examples | goes into |
|---|---|---|
| conservation laws (differential, always true) | Poisson, electron continuity, hole continuity | **PDE**, one per unknown |
| constitutive relations (local, describe the material) | n(ψ, φn) statistics, current expressions, SRH/Auger, mobility models | **Variables**, used inside c and f |
| results (follow from the above under assumptions) | Vbi, depletion width, Shockley equation | **outputs**, compared after solving |

Shockley in particular must NOT go in. It's what the model is supposed to reproduce on its own. If you put it in, you can't validate anything, and you lose the interesting part: seeing *where* it stops being true.

Fermi-Dirac isn't an extra equation either, it just changes the expression for n and p.

## why every PDE has the same shape

Almost every physics equation is a conservation law: in a small volume, change = what flows in - what flows out + what's produced inside. Divergence theorem, steady state:

$$
\nabla \cdot \boldsymbol{\Gamma} = f
$$

Γ is a flux, f a source. If the flux goes down the gradient (Fick, Fourier, Ohm), Γ = -c ∇u, and you get the simplest case of COMSOL's general form:

$$
\nabla \cdot (-c\,\nabla u - \alpha u + \gamma) + \beta \cdot \nabla u + a\,u = f
\;\;\longrightarrow\;\;
\nabla \cdot (-c\,\nabla u) = f
$$

All three diode equations fit this with only **c** and **f** non-zero. (Also why the default boundary condition is *Zero Flux*: the weak form leaves the outgoing flux on the boundary.)

## unknowns: potentials, not concentrations

n and p span 30+ orders of magnitude across the device. Solving for them directly is asking Newton to diverge. Standard trick: solve for three potentials, all in volts.

| unknown | meaning |
|---|---|
| `psi` | electrostatic potential ψ |
| `phin` | electron quasi-Fermi potential φn (E_Fn = -q φn) |
| `phip` | hole quasi-Fermi potential φp (E_Fp = -q φp) |

Concentrations become variables (Boltzmann):

$$
n = n_i \exp\!\left(\frac{\psi - \phi_n}{V_T}\right)
\qquad
p = n_i \exp\!\left(\frac{\phi_p - \psi}{V_T}\right)
$$

At equilibrium φn = φp (one Fermi level) and you get np = ni² back.

## 1. Poisson

Gauss with D = εE, E = -∇ψ:

$$
-\frac{d}{dx}\left(\varepsilon_s \frac{d\psi}{dx}\right) = q\,(p - n + N_{dop})
$$

Charge = free carriers + ionized dopants (assuming full ionization). Matching with -∇·(c∇u) = f: u = ψ, c = ε, f = ρ.

```text
interface c   (psi)
  c   = eps_s
  f   = e_const*(p-n+Ndop)
  everything else = 0
```

Units check: F/m · V/m² = C/m³ on the left, q · m⁻³ = C/m³ on the right. Yes, also in 1D (see [01](01-setup.md#what-1d-means-here)).

It's nonlinear (n, p depend on ψ). Nothing to do about that, COMSOL builds the Jacobian itself and runs Newton.

## 2. electron continuity

Current first. Drift + diffusion collapse into one term if you write it with the quasi-Fermi potential:

$$
J_n = q\mu_n n E + qD_n \frac{dn}{dx} = -q\,\mu_n\, n\, \frac{d\phi_n}{dx}
$$

??? note "derivation (click)"
    Differentiate n = nᵢ exp((ψ-φn)/Vt):

    dn/dx = (n/Vt)(dψ/dx - dφn/dx)

    plug into drift + diffusion with E = -dψ/dx and Einstein D = μ Vt:

    Jn = -qμn n dψ/dx + qμn n (dψ/dx - dφn/dx) = -qμn n dφn/dx

    The ψ terms cancel. Drift and diffusion are the same force: the slope of the electrochemical potential.

Conservation, steady state, R = net recombination rate:

$$
\frac{1}{q}\frac{dJ_n}{dx} = R
\;\;\Rightarrow\;\;
-\frac{d}{dx}\left(\mu_n\, n\, \frac{d\phi_n}{dx}\right) = R
$$

I divide both sides by ni so the numbers aren't ~1e25. Same solution, solver is happier:

```text
interface c2  (phin)
  c   = mun*n/ni
  f   = R/ni
```

## 3. hole continuity

Same thing, with J_p = -q μp p dφp/dx and (1/q) dJp/dx = -R:

$$
-\frac{d}{dx}\left(\mu_p\, p\, \frac{d\phi_p}{dx}\right) = -R
$$

```text
interface c3  (phip)
  c   = mup*p/ni
  f   = -R/ni
```

## SRH recombination

Trap at midgap:

$$
R = \frac{np - n_i^2}{\tau_p (n + n_i) + \tau_n (p + n_i)}
$$

R > 0 when np > ni² (forward bias, recombination), R < 0 in the depletion region under reverse bias (generation), R = 0 at equilibrium.

## model wizard settings

| interface | dep. variable | variable unit | source unit |
|---|---|---|---|
| Coefficient Form PDE (c) | psi | V | C/m^3 |
| Coefficient Form PDE (c2) | phin | V | 1/s |
| Coefficient Form PDE (c3) | phip | V | 1/s |

The source unit for c2/c3 is 1/s only because of the /ni scaling. Without it, 1/(m^3*s).

![model wizard](../img/screens/model_wizard.jpg)

Three Coefficient Form PDEs added in the wizard. The defaults are *Dimensionless* and *m^-2*: change both with the pencil icon, otherwise you get unit warnings everywhere (and lose an easy way to catch typos). Also rename `u`, `u2`, `u3`, mine came out as `psi3` the first time.

!!! warning "defaults"
    COMSOL sets **f = 1** by default, plus nonzero mass/damping coefficients. Change f, zero the rest. In stationary studies ea and da are ignored, but I zeroed them anyway.

## when would you need more equations

Only if a new effect needs a new unknown field:

| effect | what changes |
|---|---|
| Auger, impact ionization | extra terms in R. no new PDE |
| Fermi-Dirac, bandgap narrowing | different n, p expressions |
| self-heating | heat equation, new unknown T(x) |
| transients | time-derivative terms in the same PDEs |
| hot carriers, very short devices | energy balance / hydrodynamic |
| quantum confinement | Schrödinger-Poisson |

For a micron-sized diode, plain drift-diffusion covers the whole roadmap, breakdown included.
