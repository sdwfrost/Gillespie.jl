
<a id='Gillespie.jl:-Stochastic-Simulation-Algorithm-in-Julia-1'></a>

# Gillespie.jl: Stochastic Simulation Algorithm in Julia

- [Gillespie.jl: Stochastic Simulation Algorithm in Julia](index.md#Gillespie.jl:-Stochastic-Simulation-Algorithm-in-Julia-1)
    - [Types](index.md#Types-1)
    - [Functions](index.md#Functions-1)


<a id='Types-1'></a>

## Types

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

## Functions

<a id='Gillespie.ssa' href='#Gillespie.ssa'>#</a>
**`Gillespie.ssa`** &mdash; *Function*.



This function performs Gillespie's stochastic simulation algorithm. It takes the following arguments:

  * **x0** : a `Vector` of `Int64`, representing the initial states of the system.
  * **F** : a `Function` or a callable type, which itself takes two arguments; x, a `Vector` of `Int64` representing the states, and parms, a `Vector` of `Float64` representing the parameters of the system.
  * **nu** : a `Matrix` of `Int64`, representing the transitions of the system, organised by row.
  * **parms** : a `Vector` of `Float64` representing the parameters of the system.
  * **tf** : the final simulation time (`Float64`)

