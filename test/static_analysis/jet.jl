# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

module TestJET

using Test
using JET
using SDiagonalizability

@testset "Static analysis with JET" begin
    rep = report_package("SDiagonalizability")
    jet_reports = JET.get_reports(rep)

    @show length(jet_reports)
    @show rep

    @test length(jet_reports) < 20
    @test_broken length(jet_reports) == 0
end

end
