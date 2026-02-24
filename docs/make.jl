using Documenter, Gillespie

makedocs(
    modules = [Gillespie],
    sitename = "Gillespie.jl",
    doctest = false,
)

deploydocs(
    repo = "github.com/sdwfrost/Gillespie.jl.git",
    push_preview = true,
)
