using SDiagonalizability
using Documenter

DocMeta.setdocmeta!(
    SDiagonalizability,
    :DocTestSetup,
    :(using SDiagonalizability);
    recursive = true,
)

makedocs(;
    modules = [SDiagonalizability],
    authors = "Luis M. B. Varona <lbvarona@mta.ca>, Nathaniel Johnston <njohnston@mta.ca>",
    sitename = "SDiagonalizability.jl",
    format = Documenter.HTML(;
        canonical = "https://GraphQuantum.github.io/SDiagonalizability.jl",
        edit_link = "main",
        assets = String[],
    ),
    pages = ["Home" => "index.md"],
)

deploydocs(; repo = "github.com/GraphQuantum/SDiagonalizability.jl", devbranch = "main")
