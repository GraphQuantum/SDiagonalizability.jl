# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
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

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:alpha)

makedocs(;
    modules=[SDiagonalizability],
    authors="Luis M. B. Varona <lm.varona@outlook.com> and Nathaniel Johnston <nathaniel.johnston@gmail.com>",
    sitename="SDiagonalizability.jl",
    format=Documenter.HTML(;
        canonical="https://GraphQuantum.github.io/SDiagonalizability.jl",
        edit_link="main",
        assets=["assets/styles.css"],
        size_threshold=1_000_000,
        size_threshold_warn=200_000,
    ),
    plugins=[bib],
    pages=[
        "Home" => "index.md",
        "Public API" => "public_api.md",
        "Private API" => "private_api.md",
    ],
)

deploydocs(; repo="github.com/GraphQuantum/SDiagonalizability.jl", devbranch="main")
