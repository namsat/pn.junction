# assumptions

Everything the current model (A1) ignores, when it stops being OK, and when I plan to fix it.

| assumption | breaks down when | fixed in |
|---|---|---|
| abrupt junction | always, in real devices (diffusion/implant) | D |
| uniform doping per region | same | D |
| two regions only, no epi / n+ substrate | series resistance, punch-through matter | D |
| 1D | at the junction edges | B, C |
| Boltzmann statistics | above ~1e19 cm⁻³ (my p+ side is right at the limit) | A4 |
| complete ionization | low T, intermediate doping | probably never |
| constant bandgap | above ~1e18 cm⁻³ (bandgap narrowing) | A4 |
| constant mobilities | always a bit; a lot in the p+ and at high field | A3/A4 |
| SRH only | Auger at high doping / high injection | A4 |
| constant lifetimes | heavily doped regions | A4 |
| midgap traps (n1 = p1 = ni) | real trap levels | probably never |
| no impact ionization | fields above ~3e5 V/cm | A5 |
| uniform T = 300 K | T sweep, self-heating | A6 |
| steady state | transients, reverse recovery | maybe later |
| ideal ohmic contacts | real contacts (resistance, Schottky) | probably never |
| drift-diffusion | sub-100 nm, hot carriers | not needed for this diode |

## assumptions only the *theory* makes

The analytic formulas I compare against assume more than the simulation does. That's the useful part: wherever the two disagree, one of these is failing.

| analytic assumption | what it looks like in the simulation |
|---|---|
| depletion approximation (sharp edges) | edges smeared over a few Debye lengths, field spike at the junction (see [05](pn-1d/05-equilibrium.md#electric-field)) |
| low injection | I-V slope changes above ~0.7 V |
| no recombination in the SCR (Shockley) | ideality ≈ 2 at low bias |
| field-free neutral regions | small field in neutral regions at high current |
| no series resistance | I-V bends over at high current |
