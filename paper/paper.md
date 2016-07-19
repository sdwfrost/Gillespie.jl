---
title: 'Gillespie.jl: Stochastic Simulation Algorithm in Julia'
tags:
  - Julia
  - stochastic simulation
  - Doob-Gillespie algorithm
authors:
 - name: Simon DW Frost
   orcid: 0000-0002-5207-9879
   affiliation: University of Cambridge
date: 18 July 2016
bibliography: paper.bib
---

# Summary

`Gillespie.jl` [@GillespieJL] is a Julia package for stochastic simulation using Gillespie's direct method (sometimes called the Doob-Gillespie algorithm) [@Doob1945;@Gillespie1977] for performing stochastic simulations, an approach widely used in many fields, including systems biology and epidemiology. It borrows the basic interface (although none of the code) from the R library `GillespieSSA` by Mario Pineda-Krch [@Pineda-Krch2008], although `Gillespie.jl` only implements the standard exact method at present, whereas `GillespieSSA` also includes other methods, such as tau-leaping, *etc.*. `Gillespie.jl` is intended to offer performance on par with hand-coded C code, while maintaining a simple but flexible interface.

# References
