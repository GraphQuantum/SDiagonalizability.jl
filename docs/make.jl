# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using SDiagonalizability
using Documenter
using DocumenterCitations

DocMeta.setdocmeta!(
    SDiagonalizability, :DocTestSetup, :(using SDiagonalizability); recursive=true
)

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(;
    modules=[SDiagonalizability],
    authors="Luis M. B. Varona <lbvarona@mta.ca>, Nathaniel Johnston <njohnston@mta.ca>",
    sitename="SDiagonalizability.jl",
    format=Documenter.HTML(;
        canonical="https://GraphQuantum.github.io/SDiagonalizability.jl",
        edit_link="main",
        assets=String[],
    ),
    plugins=[bib],
    pages=["Home" => "index.md"],
)

deploydocs(; repo="github.com/GraphQuantum/SDiagonalizability.jl", devbranch="main")
