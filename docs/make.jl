using Documenter, Gillespie

makedocs()

deploydocs(
    deps = Deps.pip("mkdocs","python-markdown-math"),
    repo = "github.com/sdwfrost/Gillespie.jl.git",
    julia = "0.4"
)
