# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    TestAqua

Static analysis with `Aqua.jl` for the `SDiagonalizability.jl` package.

`Aqua.jl` offers a general static analyzer, checking for method ambiguities, undefined
`export`s, unbound type parameters, stale dependencies, type piracies, precompilation
issues, and more.
"""
module TestAqua

using SDiagonalizability
using Test
using Aqua

@testset "Static analysis with Aqua" begin
    @test Test.detect_ambiguities(SDiagonalizability) == Tuple{Method,Method}[]
    Aqua.test_all(
        SDiagonalizability;
        piracies=(; treat_as_own=[SDiagonalizability.LinearAlgebra.rank]),
        persistent_tasks=false, # Account for our manual definition of `LinearAlgebra.rank`
    )
    #= We manually define the `LinearAlgebra.rank(::QRPivoted)` method, as it is not
    available in Julua 1.11 and earlier. =#
    @test length(Aqua.Piracy.hunt(SDiagonalizability)) <= 2
end

end
