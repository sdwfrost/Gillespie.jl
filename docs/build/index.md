
<a id='Gillespie.jl:-Stochastic-Simulation-Algorithm-in-Julia-1'></a>

# Gillespie.jl: Stochastic Simulation Algorithm in Julia

- [Gillespie.jl: Stochastic Simulation Algorithm in Julia](index.md#Gillespie.jl:-Stochastic-Simulation-Algorithm-in-Julia-1)
    - [Introduction](index.md#Introduction-1)
    - [Installation](index.md#Installation-1)
    - [Example](index.md#Example-1)
    - [Application programming interface](index.md#Application-programming-interface-1)


<a id='Introduction-1'></a>

## Introduction


`Gillespie.jl` provides an implementation of [Gillespie's direct method](http://en.wikipedia.org/wiki/Gillespie_algorithm) for performing stochastic simulations, which are widely used in many fields, including systems biology and epidemiology. It borrows the basic interface (although none of the code) from the R library [`GillespieSSA`](http://www.jstatsoft.org/v25/i12/paper) by Mario Pineda-Krch, although `Gillespie.jl` only implements the standard exact method at present, whereas `GillespieSSA` also includes tau-leaping, *etc.*.


<a id='Installation-1'></a>

## Installation


`Gillespie.jl` can be installed from the Julia read-eval-print-loop (REPL) as follows.


```julia
Pkg.add("Gillespie")
```


<a id='Example-1'></a>

## Example


Let's take the 'standard' [susceptible-infected-recovered (SIR) model]((https://en.wikipedia.org/wiki/Compartmental_models_in_epidemiology#The_SIR_model_without_vital_dynamics)), commonly used in epidemiology. A deterministic version of this model is as follows.


\begin{align} \frac{dS(t)}{dt}  & = -\beta S(t) I(t) \cr \frac{dI(t)}{dt}  & = \beta S(t) I(t)- \gamma I(t) \cr \frac{dR(t)}{dt}  & = \gamma I(t) \end{align}


Let's consider a stochastic version of the SIR model.


\begin{align} {\rm Transition} & \quad {\rm Rate} \cr S  \rightarrow S-1,\; I \rightarrow I+1 & \quad \beta S(t) I(t) \cr I  \rightarrow I-1,\; R \rightarrow R+1 & \quad \gamma I(t) \end{align}


We first need to load the library.


```julia
using Gillespie;
```


We next need to define a function that given state variables `x` (type: `Array{Int64,1}`) and a vector of parameters (type: `Vector{Float64}`), returns a vector of rates of length `k` for different types of transitions. For this example, there are two transition functions, corresponding to infection and recovery.


```julia
function F(x,parms)
  (S,I,R) = x
  (beta,gamma) = parms
  infection = beta*S*I
  recovery = gamma*I
  [infection,recovery]
end;
```

```
F (generic function with 1 method)
```


We define the states of the system - a `Vector{Int64}` of length `n` with the number of susceptible, infected, and recovered individuals.


```julia
x0 = [9999,1,0];
```

```
3-element Array{Int64,1}:
 9999
    1
    0
```


To define the transitions, we define an `Array{Int64,k,n}` that denotes the changes to each of the `n` state variables for each of the `k` transitions. Infection results in a loss of 1 susceptible and a gain of one infected individual, while recovery is associated with a loss of one infected and a gain of one recovered.


```julia
nu = [[-1 1 0];[0 -1 1]];
```

```
2x3 Array{Int64,2}:
 -1   1  0
  0  -1  1
```


Finally, we define the parameter values (in the order required by the function `F`), and the time we want the simulation to finish. The simulation will finish early if the propensity rates are zero.


```julia
parms = [0.1/10000.0,0.05]
tf = 1000.0;
```

```
1000.0
```


Given the above, the simulation can be run using the function `ssa`. It's usually a good idea to set a random number seed prior to simulation first.


```julia
using Random
Random.seed!(1236)
result = ssa(x0,F,nu,parms,tf);
```

```
Gillespie.SSAResult([0.0,11.7409,12.1425,17.647,18.0987,21.2081,22.2022,22.6326,26.1136,27.201  …  413.431,413.962,414.973,420.548,422.656,424.17,424.486,430.513,433.936,482.663],15704x3 Array{Int64,2}:
 9998  2     0
 9998  2     0
 9997  3     0
 9996  4     0
 9995  5     0
 9994  6     0
 9994  5     1
 9994  4     2
 9994  3     3
 9993  4     3
    ⋮
 2149  6  7845
 2149  5  7846
 2149  4  7847
 2148  5  7847
 2148  4  7848
 2148  3  7849
 2148  2  7850
 2148  1  7851
 2148  0  7852,Gillespie.SSAStats("zeroprop",15703),Gillespie.SSAArgs([9999,1,0],ex-1.F,2x3 Array{Int64,2}:
 -1   1  0
  0  -1  1,[1.0e-5,0.05],1000.0))
```


This will return an object of type `SSAresult`. This can be converted to a `DataFrame` using the function `ssa_data`.


```julia
data = ssa_data(result);
```

```
15704×4 DataFrames.DataFrame
│ Row   │ time    │ x1   │ x2 │ x3   │
├───────┼─────────┼──────┼────┼──────┤
│ 1     │ 0.0     │ 9998 │ 2  │ 0    │
│ 2     │ 11.7409 │ 9998 │ 2  │ 0    │
│ 3     │ 12.1425 │ 9997 │ 3  │ 0    │
│ 4     │ 17.647  │ 9996 │ 4  │ 0    │
│ 5     │ 18.0987 │ 9995 │ 5  │ 0    │
│ 6     │ 21.2081 │ 9994 │ 6  │ 0    │
│ 7     │ 22.2022 │ 9994 │ 5  │ 1    │
│ 8     │ 22.6326 │ 9994 │ 4  │ 2    │
⋮
│ 15696 │ 413.962 │ 2149 │ 6  │ 7845 │
│ 15697 │ 414.973 │ 2149 │ 5  │ 7846 │
│ 15698 │ 420.548 │ 2149 │ 4  │ 7847 │
│ 15699 │ 422.656 │ 2148 │ 5  │ 7847 │
│ 15700 │ 424.17  │ 2148 │ 4  │ 7848 │
│ 15701 │ 424.486 │ 2148 │ 3  │ 7849 │
│ 15702 │ 430.513 │ 2148 │ 2  │ 7850 │
│ 15703 │ 433.936 │ 2148 │ 1  │ 7851 │
│ 15704 │ 482.663 │ 2148 │ 0  │ 7852 │
```


This makes it straightforward to plot e.g. using `Gadfly`.


```julia
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
```


![](plot.svg)


<a id='Application-programming-interface-1'></a>

## Application programming interface


<a id='Types-1'></a>

### Types

<a id='Gillespie.SSAStats' href='#Gillespie.SSAStats'>#</a>
**`Gillespie.SSAStats`** &mdash; *Type*.



A type storing the status at the end of a call to `ssa`.

<a id='Gillespie.SSAArgs' href='#Gillespie.SSAArgs'>#</a>
**`Gillespie.SSAArgs`** &mdash; *Type*.



A type storing the call to `ssa`.

<a id='Gillespie.SSAResult' href='#Gillespie.SSAResult'>#</a>
**`Gillespie.SSAResult`** &mdash; *Type*.



This type stores the output of `ssa`, and comprises of:

  * **time** : a `Vector` of `Float64`, containing the times of simulated events.
  * **data** : a `Matrix` of `Int64`, containing the simulated states.
  * **stats** : an instance of `SSAStats`.
  * **args** : arguments passed to `ssa`.


<a id='Functions-1'></a>

### Functions

<a id='Gillespie.ssa' href='#Gillespie.ssa'>#</a>
**`Gillespie.ssa`** &mdash; *Function*.



This function performs Gillespie's stochastic simulation algorithm. It takes the following arguments:

  * **x0** : a `Vector` of `Int64`, representing the initial states of the system.
  * **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
  * **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
  * **parms** : a `Vector` of `Float64` representing the parameters of the system.
  * **tf** : the final simulation time (`Float64`)

<a id='Gillespie.ssa_data' href='#Gillespie.ssa_data'>#</a>
**`Gillespie.ssa_data`** &mdash; *Function*.



This takes a single argument of type `SSAResult` and returns a `DataFrame`.


```
pfsample
```

