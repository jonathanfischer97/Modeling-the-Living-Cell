using Random 
using Plots 
using Catalyst
using DifferentialEquations
using Latexify
using Graphviz_jll
# using DiffEqBiological 

 #main model. Gene network from HW7 with the condensation mechanism added 
osc_rn = @reaction_network osc begin
    kf1, geneA --> geneA + mRNA_A  
    (kf2,kb2), geneA + A <--> geneA_bound
    kf3, geneA_bound --> geneA_bound + mRNA_A 
    kf4, geneR --> geneR + mRNA_R 
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
    kf15, A --> Ag #condensation mechanism. A aggregates to form Ag
    kf16, Ag + X --> A + X #phantom species X scales the rate of disassembly by the V/A ratio 
    kf17, Ag + X --> A + Xseq #another phantom reaction to preserve mass conservation 
    kf18, Xseq --> X
    kf19, geneR + Ag --> geneR_Ag #Superenhancer aggregate binds to geneR to form geneR_Ag
    kf20, geneR_Ag --> geneR + Ag 
    kf21, geneR_Ag --> geneR_Ag + mRNA_R #Enhancer stimulates production of mRNA_R
end kf1 kf2 kb2 kf3 kf4 kf5 kb5 kf6 kf7 kf8 kf9 kf10 kf11 kf12 kf13 kf14 kf15 kf16 kf17 kf18 kf19 kf20 kf21  
osys = convert(ODESystem, osc_rn)
copy_to_clipboard(false)
println(latexify(osc_rn; env=:arrow))

 #this model was unused in the presentation, but was a prototype for the later models 
drop_rn = @reaction_network drop begin
    #D is produced in the cytosol through a bimolecular reaction of A and B, and is consumed in the droplet when A binds to an enzyme E to catalyze the degradation of D. Oscillations in the size of the droplet are produced because the size of the droplet determines the rate at which A and B can diffuse into the droplet to activate the enzyme. A and B can diffuse back out of the droplet and into the cytosol according to the law of mass action.
    #The rate of the diffusion reactions is dependent on the size of the droplet, which is determined by the concentration of D in the droplet. The concentration of D in the droplet is determined by the rate of the bimolecular reaction of A and B in the cytosol, which is dependent on the concentration of A and B in the cytosol. The concentration of A and B in the cytosol is determined by the rate of the diffusion reactions, which is dependent on the size of the droplet. Thus, the size of the droplet oscillates between two values, which is the basis for the oscillations in the concentration of A and B in the cytosol.
    #The diffusive rate constants are scaled by Vdrop/Vref, or the ratio of the droplet volume (variable) to a reference droplet size (fixed parameter)
    (ka1,kb1), A_cyt + B_cyt <--> D_cyt # A and B bind to form D in the cytosol 
    (ka2*Vdrop/Vref, kb2), D_cyt <--> D_drop # D diffuses and moves to the drop. Forward rate is scaled by the size of the droplet, which is a function of its total concentration. 
    (ka3, kb3), D_drop + EA <--> D_dropEA # D and the activated enzyme EA bind to form intermediate complex 
    kcat3, D_dropEA --> A_drop + B_drop + EA # EA catalyzes the hydrolysis of D to form A and B in the droplet 
    (ka4*Vdrop/Vref, kb4), A_drop <--> A_cyt # A diffuses and moves to the cytoplasm
    (ka5*Vdrop/Vref, kb5), B_drop <--> B_cyt # B diffuses and moves to the cytoplasm
    (ka6, kb6), E + A_drop <--> EA # E binds to A to form EA in the droplet 
    (kd*A_drop*B_drop*D_drop, kd*A_cyt*B_cyt*D_cyt), 0 <-->  Vdrop # Reaction to relate total droplet concentration to droplet volume with proportionality constant kd 
end ka1 kb1 ka2 kb2 ka3 kb3 kcat3 ka4 kb4 ka5 kb5 ka6 kb6 kd Vref 
osys2 = convert(ODESystem, drop_rn)
println(latexify(drop_rn; env=:arrow))


#function to calculate surface area to volume ratio given a volume
function AV_ratio(volume)
    return (4*pi*(volume/3)^(2/3))/volume
end

#function to calculate volume of a sphere given a copy number and estimated bond length 
function get_volume(copies, bondlength=1)
    radius = bondlength/2 
    return copies*4/3*pi*radius^3
end

"AV = AV_ratio(get_volume(Ag,bondlength))" #This is how the A/V ratio relates to copy numbers of Ag. The bondlength is the length of the bond between Ag molecules. The default value is 1 um.

#This model is the idealized version taking advantage of condensate positive feedback driving kinase mediated negative feedback. Didn't work for optimization. 
drop_rn2 = @reaction_network drop2 begin 
    (ka1, kb1), Ap + B <--> ApB #only phosphorylated A can dimerize 
    (ka2, kb2*AV), ApB <--> ApBg #Ag stands for aggregate of A dimers 
    (ka3, kb3*AV), K + X <--> Kseq + X #Kinase is sequestered in droplet 
    (ka4, kb4), K + A <--> KA #Kinase binds to A 
    kcat4, KA --> K + Ap #Kinase catalyzes the phosphorylation of A
    (ka5, kb5), P + Ap <--> PAp #Phosphatase binds to Ap
    kcat5, PAp --> P + A #Phosphatase catalyzes the dephosphorylation of Ap
end ka1 kb1 ka2 kb2 ka3 kb3 ka4 kb4 kcat4 ka5 kb5 kcat5 AV 
osys3 = convert(ODESystem, drop_rn2)
println(latexify(drop_rn2; env=:arrow))


#Here is the final model before the main model, using phantom species and reactions to temper the ultrasensitivity of the system.
drop_rn3 = @reaction_network drop3 begin 
    (ka1, kb1), Ap + B <--> ABdimer #only phosphorylated A can dimerize with B
    ka2, ABdimer --> ABg #ABg stands for aggregate of AB heterodimer aggregation 
    kb2, ABg + X --> ABdimer + X #Phantom species X scales disassembly of ABg
    (ka3, kb3), K + ABg <--> Kseq + ABg #Kinase is sequestered in droplet, driving long negative feedback 
    (ka4, kb4), K + A <--> KA #Kinase binds to A 
    kcat4, KA --> K + Ap #Kinase catalyzes the phosphorylation of A
    (ka5, kb5), P + Ap <--> PAp #Phosphatase binds to Ap
    kcat5, PAp --> P + A #Phosphatase catalyzes the dephosphorylation of Ap
    ka6, ABg + X --> Xseq  + ABg #Mass conservation phantom reaction 
    kb6, Xseq --> X
end ka1 kb1 ka2 kb2 ka3 kb3 ka4 kb4 kcat4 ka5 kb5 kcat5 ka6 kb6 
osys3 = convert(ODESystem, drop_rn3)
println(latexify(drop_rn3; env=:arrow))


# Print out the ODEs of the MAIN MODEL without the tildes and parentheses for use in genetic algorithm in Python 
for (i,ode) in enumerate(osys.eqs)
    ode_str = string(ode)
    ode_str = replace(ode_str, "~" => "=")
    ode_str = replace(ode_str, "(t)" => "")
    #replace Differential(var) with dvar 
    ode_str = replace(ode_str, "Differential(" => "d") 
    #remove other paranthesis on the left hand side of the equation
    ode_str = replace(ode_str, ") =" => " =")
    println(ode_str)
end


#optimized parameter values from genetic algorithm
p = [18.74707803532138, 13.358801863610093, 200.0, 200.0, 41.47638690667699, 92.40925991373919, 62.51404128179669, 27.778548458556934, 164.56149941287393, 75.75725724994625, 133.8004596066915, 3.3942171322224852, 37.56958435766569, 165.1953424131399, 2.025279604137535, 7.846378991546339, 49.28084335290114, 14.002254763294918, 45.681547362024055, 47.456140095892735, 76.45946322667373, 12.943552212252534, 173.99428155159933]
pnames = [:kf1, :kf2, :kb2, :kf3, :kf4, :kf5, :kb5, :kf6, :kf7, :kf8, :kf9, :kf10, :kf11, :kf12, :kf13, :kf14, :kf15, :kf16, :kf17, :kf18, :kf19, :kf20, :kf21]
pmap = [x => y for (x,y) in zip(pnames, p)]
pmap = symmap_to_varmap(osc_rn, pmap)

#Paired list of initial conditions 
umap = [:geneA => 1., :mRNA_A => 0., :A => 0., :geneA_bound => 0., :geneR => 1., :mRNA_R => 0., :geneR_bound => 0., :R => 0., :C => 0., :Ag => 0., :X => 100., :Xseq => 0., :geneR_Ag => 0.]
umap = symmap_to_varmap(osc_rn, umap)




tspan = (0.,100.)
oprob = ODEProblem(osc_rn, umap, tspan, pmap)
osol = solve(oprob)



#Plot time series of oscillations for A, R, C, and Ag
plot(osol.t, osol[3,:], color = :blue, label = "A")
plot!(osol.t, osol[8,:], color = :green, label = "R")
plot!(osol.t, osol[9,:], color = :pink, label = "C")
plot!(osol.t, osol[10,:], color = :gray, line = :dash, label = "Ag")
xaxis!("Time (s)")
yaxis!("Copies")
savefig("./osc.png")

#Plot time series for the volume of the droplet as function of Ag copies 
plot(osol.t, get_volume(osol[10,:]), color = :green, label = "Condesate Volume")
xaxis!("Time (s)")
yaxis!("Volume (nm^3)")
savefig("./osc_volume.png")



# Graph reaction network osc_rn
g = Graph(osc_rn)
savegraph(g, "./drop_rn_graph.png")




#plot phase diagram of Ag vs C with axes labeled as copy numbers 
plot(get_volume(osol[10,:]),osol[9,:], osol[3,:], color = :blue, label = "C vs Volume vs A")
savefig("./phase.png")
