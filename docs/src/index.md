# Gillespie.jl: Stochastic Simulation Algorithm in Julia

```@contents
```

## Introduction

`Gillespie.jl` provides an implementation of [Gillespie's direct method](http://en.wikipedia.org/wiki/Gillespie_algorithm) for performing stochastic simulations, which are widely used in many fields, including systems biology and epidemiology. It borrows the basic interface (although none of the code) from the R library [`GillespieSSA`](http://www.jstatsoft.org/v25/i12/paper) by Mario Pineda-Krch, although `Gillespie.jl` only implements the standard exact method at present, whereas `GillespieSSA` also includes tau-leaping, *etc.*.

## Installation

`Gillespie.jl` can be installed from the Julia read-eval-print-loop (REPL) as follows.

```julia
Pkg.add("Gillespie")
```

## Example

Let's take the 'standard' [susceptible-infected-recovered (SIR) model]((https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model_without_vital_dynamics)), commonly used in epidemiology. A deterministic version of this model is as follows.

\begin{align}
\frac{dS(t)}{dt}  & = -\beta S(t) I(t) \cr
\frac{dI(t)}{dt}  & = \beta S(t) I(t)- \gamma I(t) \cr
\frac{dR(t)}{dt}  & = \gamma I(t)
\end{align}

Let's consider a stochastic version of the SIR model.

\begin{align}
{\rm Transition} & \quad {\rm Rate} \cr
S  \rightarrow S-1,\; I \rightarrow I+1 & \quad \beta S(t) I(t) \cr
I  \rightarrow I-1,\; R \rightarrow R+1 & \quad \gamma I(t)
\end{align}

We first need to load the library.

```@example 1
using Gillespie
```

We next need to define a function that given state variables `x` (type: `Array{Int64,1}`) and a vector of parameters (type: `Vector{Float64}`), returns a vector of rates of length `k` for different types of transitions. For this example, there are two transition functions, corresponding to infection and recovery.

```@example 1
function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end
```

We define the states of the system - a `Vector{Int64}` of length `n` with the number of susceptible, infected, and recovered individuals.

```@example 1
x0 = [9999,1,0]
```

To define the transitions, we define an `Array{Int64,k,n}` that denotes the changes to each of the `n` state variables for each of the `k` transitions. Infection results in a loss of 1 susceptible and a gain of one infected individual, while recovery is associated with a loss of one infected and a gain of one recovered.

```@example 1
nu = [[-1 1 0];[0 -1 1]]
```

Finally, we define the parameter values (in the order required by the function `F`), and the time we want the simulation to finish. The simulation will finish early if the propensity rates are zero.

```@example 1
parms = [0.1/10000.0,0.05]
tf = 1000.0
```

Given the above, the simulation can be run using the function `ssa`. It's usually a good idea to set a random number seed prior to simulation first.

```@example 1
srand(1234)
result = ssa(x0,F,nu,parms,tf)
```

This will return an object of type `SSAresult`. This can be converted to a `DataFrame` using the function `ssa_data`.

```@example 1
data = ssa_data(result)
```

This makes it straightforward to plot e.g. using `Gadfly`.

```@example 1
using Gadfly
plot(data,
  layer(x="time",y="x1",Geom.step,Theme(default_color=colorant"red")),
  layer(x="time",y="x2",Geom.step,Theme(default_color=colorant"blue")),
  layer(x="time",y="x3",Geom.step,Theme(default_color=colorant"green")),
  Guide.xlabel("Time"),
  Guide.ylabel("Number"),
  Guide.manual_color_key("Population",
                            ["S", "I", "R"],
                            ["red", "blue", "green"]),
  Guide.title("SIR epidemiological model"))
draw(SVG("plot.svg",6inch,4inch),ans); nothing # hide
```

![](plot.svg)

## Application programming interface

### Types

```@docs
SSAStats
```

```@docs
SSAArgs
```

```@docs
SSAResult
```

### Functions

```@docs
ssa
```

```@docs
ssa_data
```
