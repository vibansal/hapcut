
HAPCUT code (designed to handle long reads) that does not use the read-haplotype graph for max-cut computations but uses a greedy heuristic to update
the haplotypes using the exact haplotype likelihoods.

Features of this code:

1. Operates on likelihood function and also models chimeric fragments
2. storage requirements for graph are significantly less for long reads 
