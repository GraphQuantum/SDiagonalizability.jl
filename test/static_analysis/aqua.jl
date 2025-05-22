# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module TestAqua

using Test
using Aqua
using SDiagonalizability

@testset "Static analysis with Aqua" begin
    @test Test.detect_ambiguities(SDiagonalizability) == Tuple{Method,Method}[]
    Aqua.test_all(SDiagonalizability)
    @test Aqua.Piracy.hunt(SDiagonalizability) == Method[]
end

end
