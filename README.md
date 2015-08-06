# Gillespie

[![Build Status](https://travis-ci.org/sdwfrost/Gillespie.jl.svg?branch=master)](https://travis-ci.org/sdwfrost/Gillespie.jl)

This is a preliminary implementation of [Gillespie's direct method](http://en.wikipedia.org/wiki/Gillespie_algorithm) for performing stochastic simulations. It borrows the basic interface (although none of the code) from the R library [`GillespieSSA`](http://www.jstatsoft.org/v25/i12/paper) by Mario Pineda-Krch, although `Gillespie.jl` only implements the standard exact method at present, whereas `GillespieSSA` also includes tau-leaping, *etc.*.
