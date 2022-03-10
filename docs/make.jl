using GSA_CME
using Documenter

DocMeta.setdocmeta!(GSA_CME, :DocTestSetup, :(using GSA_CME); recursive=true)

makedocs(;
    modules=[GSA_CME],
    authors="Aniket Jivani",
    repo="https://gitlab.umich.edu/ajivani/GSA_CME.jl/blob/{commit}{path}#{line}",
    sitename="GSA_CME.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://ajivani.gitlab.io/GSA_CME.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
