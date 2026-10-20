# devlog

## 2026-10-20

- decided on the plan: pn first (1D → 2D → 3D), BJT and MOS later
- no Semiconductor Module → writing drift-diffusion by hand with 3 Coefficient Form PDEs
- went back and forth on the structure: started with p+ / n-epi / n+ substrate, dropped it. Too much at once. Simple abrupt p+/n first, realistic structure in phase D
- model wizard: 3 interfaces, fixed units (V, C/m^3, 1/s)
- geometry, parameters, variables, BCs
- errors: `Undefined variable` x3 (Variables node on a point, missing `ni`, missing `eps_s`)
- first ψ plot wrong at the cathode (wrong Dirichlet expression). Fixed
- **equilibrium works**: Vbi ≈ 0.89 V, matches
- wrote `pn1d_reference.py` to have independent numbers. Found an `np.isclose` tolerance bug in my own script (1e-8 absolute = 10 nm here...). After fixing it there's still a field spike at the junction, turns out it's physical (hole spillover, abrupt profile)

next:

- [ ] zoom plots in COMSOL, measure W and Emax
- [ ] export COMSOL line plots to CSV, overlay with python
- [ ] step 2: voltage sweep, first I-V
