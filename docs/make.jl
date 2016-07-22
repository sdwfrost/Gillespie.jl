using Documenter, Gillespie

makedocs(modules=[Gillespie],
    doctest = false)

deploydocs(
    deps = Deps.pip("mkdocs","python-markdown-math"),
    repo = "github.com/sdwfrost/Gillespie.jl.git",
    julia = "0.4.6",
    osname = "linux"
)
