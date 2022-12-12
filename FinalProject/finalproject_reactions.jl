using Random 
using Plots 
using Catalyst
using DifferentialEquations
using Combinatorics
using Latexify

osc_rn = @reaction_network osc begin
    kf1, geneA --> geneA + mRNA_A 
    (kf2,kb2), geneA + A <--> geneA_bound
    kf3, geneA_bound --> geneA_bound + mRNA_A 
    kf4, geneR --> mRNA_R + mRNA_R 
    (kf5,kb5), geneR + A <--> geneR_bound
    kf6, geneR_bound --> geneR_bound + mRNA_R
    kf7, mRNA_A --> mRNA_A + A
    kf8, mRNA_R --> mRNA_R + R
    kf9, A + R --> C 
    kf10, C --> R  
    kf11, A --> ∅
    kf12, R --> ∅
    kf13, mRNA_A --> ∅
    kf14, mRNA_R --> ∅
end kf1 kf2 kb2 kf3 kf4 kf5 kb5 kf6 kf7 kf8 kf9 kf10 kf11 kf12 kf13 kf14   

osys = convert(ODESystem, osc_rn)

latexify(osys.eqs)

osys.eqs[2]

for eq in osys.eqs 
    println(eq)
end

typeof(osys.eqs[1])