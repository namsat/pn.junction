# pn.junction

Building semiconductor devices in COMSOL **without the Semiconductor Module**. Poisson + the two continuity equations, typed in by hand as *Coefficient Form PDEs*, then checked against textbook formulas and a small Python solver I wrote on the side.

Started with a pn diode. BJT and MOSFET come later (maybe).

> **status:** work in progress. 1D equilibrium works. Current/I-V is next. Expect broken stuff and TODOs.

## why no Semiconductor Module

Mostly because I don't have it. But it turned out to be the better way to learn this:

- you see every equation that gets solved, nothing is hidden behind a checkbox
- every assumption (Boltzmann statistics, constant mobility, SRH only...) is something *you* wrote, so you can also remove it later
- the same setup works in any PDE solver, not just COMSOL

## where to start

If you want to reproduce the model, go in order:

1. [setup + parameters](pn-1d/01-setup.md) - geometry, doping, the parameter table
2. [equations](pn-1d/02-equations.md) - drift-diffusion, and how it maps to COMSOL's coefficient form
3. [boundary conditions](pn-1d/03-boundary-conditions.md) - where the contact values actually come from
4. [mesh + solver](pn-1d/04-mesh-solver.md) - mesh grading, solver steps, the errors I ran into
5. [equilibrium results](pn-1d/05-equilibrium.md) - plots, numbers, comparison with theory

The [assumptions](assumptions.md) page keeps track of everything I'm ignoring right now and when I plan to stop ignoring it. The plan for the whole thing is in [roadmap](roadmap.md).

## repo layout

```text
pn.junction/
├── analytics/pn1d_reference.py   # analytic numbers + independent Poisson solve
├── comsol/                       # .mph files (see comsol/README.md)
├── data/reference/               # expected values (json), python solution (csv)
├── docs/                         # this site (mkdocs)
└── mkdocs.yml
```

## running the reference script

```bash
pip install numpy matplotlib
python analytics/pn1d_reference.py
```

It regenerates every plot on this site and dumps the expected values to `data/reference/expected_values.json`.

## running the site locally

```bash
pip install mkdocs-material
mkdocs serve
```
