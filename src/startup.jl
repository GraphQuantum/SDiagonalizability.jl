# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

@setup_workload begin
    g = complete_graph(6)
    L_g = laplacian_matrix(g)

    h = cartesian_product(complete_graph(2), complete_graph(4))
    L_h = laplacian_matrix(h)

    @compile_workload begin
        for network in [g, L_g, h, L_h]
            for S in [(-1, 0, 1), (-1, 1)]
                s_bandwidth(network, S)

                for k in [1, 2, 5]
                    has_s_bandwidth_at_most_k(network, S, k)
                end

                is_s_diagonalizable(network, S)
            end
        end
    end
end
