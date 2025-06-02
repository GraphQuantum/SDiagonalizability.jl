eigvals_sorted = collect(Iterators.flatmap(pair -> fill(pair...), eigval_counts))
