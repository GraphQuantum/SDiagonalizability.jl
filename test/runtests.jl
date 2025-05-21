# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

using Test
using SDiagonalizability

for analyzer in readlines(joinpath(@__DIR__, "staticanalyzers"))
    @info "Running static analysis with $analyzer"
    include("static_analysis/$(lowercase(analyzer)).jl")
    println()
end

for file in readlines(joinpath(@__DIR__, "testgroups"))
    @info "Testing $file"
    include("$file.jl")
    println()
end

@testset "Docstrings" begin
    # TODO: Change from `@test_broken` to `@test` when the issue is fixed
    @test_broken isempty(Docs.undocumented_names(SDiagonalizability))
end
