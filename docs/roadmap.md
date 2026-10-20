# roadmap

Rule: don't move to the next step until the current one matches theory (or I understand *why* it doesn't).

Device for all of phase A: abrupt p+/n, NA = 1e19, ND = 1e16 cm⁻³.

## A - 1D

| step | what | check against | status |
|---|---|---|---|
| A0 | parameters, geometry, doping, contacts, mesh | doping profile looks right | done |
| A1 | equilibrium, Poisson only | Vbi, W, Emax | done |
| A2 | I-V with long lifetimes | Shockley J0 | next |
| A3 | realistic τ, sweep τ | n ≈ 2 → 1, SRH peak in SCR, reverse current ~ √V | todo |
| A4 | Auger, bandgap narrowing, Fermi-Dirac (one at a time) | how much each one moves the I-V | todo |
| A5 | breakdown (impact ionization) | BV ≈ 50-60 V (Sze) | todo |
| A6 | T = 250-400 K | ≈ -2 mV/K forward, BV goes up with T | todo |
| A7 | small-signal C-V | ND and Vbi from 1/C² | todo |

## B - 2D axisymmetric

Circular diode, p+ diffused into n with some lateral spread. Only for edge effects:

- [ ] B1 cut through the center = 1D result (sanity check)
- [ ] B2 field crowding at the curved junction edge
- [ ] B3 breakdown at the edge (should be lower than 1D, compare with cylindrical junction formula, Baliga)
- [ ] B4 area vs perimeter current/capacitance
- [ ] B5 floating guard ring (optional)

## C - 3D, quarter model

- [ ] C1 check against 2D along the edge
- [ ] C2 breakdown at the corner (spherical junction, should be even lower)
- [ ] C3 rounded corner radius sweep (optional)

## D - realistic diode

- [ ] gaussian p+ profile instead of a step
- [ ] n epi on n+ substrate, series resistance, punch-through vs epi thickness
- [ ] design for a spec (BV = 100 V), size it by hand, verify

## E - compact model

- [ ] extract Is, n, Rs, Cj0, Vj, M, BV
- [ ] drop them into ngspice / Simscape, simple rectifier
