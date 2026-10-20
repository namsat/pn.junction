# pn.junction

pn diode (and later BJT / MOSFET) in COMSOL **without the Semiconductor Module**: drift-diffusion written by hand as three Coefficient Form PDEs, checked against textbook formulas and a small independent Python solver.

Site / tutorial: https://namsat.github.io/pn.junction/

**Status:** WIP. 1D equilibrium done, I-V next.

```
analytics/   reference numbers + independent Poisson solver (python)
comsol/      .mph models
data/        expected values, python solution, COMSOL exports
docs/        the site (mkdocs-material)
```

```bash
pip install -r requirements.txt
python analytics/pn1d_reference.py   # regenerates plots + data/reference/
mkdocs serve                         # site on localhost:8000
```
