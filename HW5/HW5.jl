using Random, Distributions
using Plots
# plotlyjs()
# using PlotlyJS
using DataFrames
using StatsBase
using DelimitedFiles
using LinearAlgebra
# using Distances
# using Symbolics
# using Latexify

function initR()
    R = zeros((5,3))
    R[:,3] .= 1,2,3,4,5
    R
end

function calc_dist(x1,x2) #calculate distance between a single pair of points along single component
    x1-x2 
    # dx - L*round(dx/L) #distance using minimum image convention 
end

function calc_LJ(r, epsilon, sigma = 1) # calculate Lennard_Jones potential for a pairwise distance
    4*epsilon*((sigma/r)^12 - (sigma/r)^6) 
end

function calc_LJforce(r, epsilon, T, v, dt, sigma = 1, k=1, m=1, eta=0.05)
    # (48*epsilon*((sigma^12/r^13)-0.5*(sigma^6/r^7)) + randn()*√(2kT*m*eta/dt) - eta*m*v, 
    (48*epsilon*((sigma^12/r^13)-0.5*(sigma^6/r^7)),
    randn()*√(2k*T*m*eta/dt),
    -(eta*m*v))
end

function calc_harmU(r,k=20,l=1)
    0.5k*(r-l)^2
end

function calc_harmforce(r,k=20,l=1)
    k*(r-l),0.,0.
end

function update_A(R, T, V, dt, m = 1) #update distance matrix D and acceleration matrix A 
    PE = 0 #initialize potential energy
    N = 5
    A = zeros((N,3,3))
    for i in 1:N-1
        for j in i+1:N
            dr = calc_dist.(R[i,:],R[j,:]) #calculate distance vector between two points
            r_mag = norm(dr) #calculate magnitude of distance vector
            if j - i > 1
                if (i,j) == (1,4) || (i,j) == (1,5) || (i,j) == (4,5)
                    epsilon = 1.
                    PE += calc_LJ(r_mag, epsilon)  #accumulate LJ potential energy
                    # f_scalars = calc_LJforce(r_mag,epsilon,T,v,dt) #scalar magnitude of force
                    A[i,:,:] += [x * dr[ii] / r_mag for x in calc_LJforce(r_mag,epsilon,T,norm(V[i]),dt), ii in 1:3] #force vector with each component
                    A[j,:,:] -= [x * dr[ii] / r_mag for x in calc_LJforce(r_mag,epsilon,T,norm(V[j]),dt), ii in 1:3] #force vector with each component
                else
                    epsilon = 2/3
                    PE += calc_LJ(r_mag, epsilon)  #accumulate LJ potential energy
                    A[i,:,:] += [x * dr[ii] / r_mag for x in calc_LJforce(r_mag,epsilon,T,norm(V[i]),dt), ii in 1:3] #force vector with each component
                    A[j,:,:] -= [x * dr[ii] / r_mag for x in calc_LJforce(r_mag,epsilon,T,norm(V[j]),dt), ii in 1:3] #force vector with each component
                end
            else 
                PE += calc_harmU(r_mag)
                # fvector = [x * dr/ r_mag for x in calc_harmforce(r_mag)]
                fvector = [x * dr[ii] / r_mag for x in calc_harmforce(r_mag), ii in 1:3]
                # return fvector
                A[i,:,:] += fvector
                A[j,:,:] -= fvector
            end
        end
    end
    return PE, A
end

function calc_position(r, v, a, dt) #update position vector for single particle
    r + (v*dt) + 0.5(a*dt^2)
end

function update_R!(R, V, A, dt) #updates position matrix R IN PLACE
    broadcast!(calc_position,R,R,V,A,dt) #vectorized operation passed R by reference 
end

function calc_velocity(v,a1,a2c,a2r,dt,eta=0.05) #velocity verlet update of velocity for single component
    (1+dt*eta/2)^-2 * (v + 0.5(a1)*dt + 0.5(a2c + a2r)*dt) #average current and new accelerations
end

function update_V!(V,A1,A2c,A2r,dt) #updates velocity matrix V IN PLACE (passed by reference)
    # A2c = copy(A2[:,:,1])
    # A2r = copy(A2[:,:,2])
    broadcast!(calc_velocity,V,V,A1,A2c,A2r,dt)
end

function init_velocity(T, N, m=1, k=1) #temperature, count, and mass
    V = √(k*T/m)*randn((N,3)) #return N x 3 matrix of initialized velocities sampled from Maxwell-Boltzmann distribution
    V[end,:] = -sum(V[1:end-1,:],dims=1) #Momentum conservation: last row equals minus sum of all other rows. Mass cancels out because same species 
    return V
end





function calc_KE(V, N, m =1)
    KE = 0.
    T = 0.
    for p in eachrow(V)
        v_squared = sum(p.^2) #get squared magnitude of velocity vector
        KE += 0.5(m*v_squared) #accumulate KE for each particle
        T += v_squared/(3N-3)
    end
    return KE, T 
end




function LD(T, dt, iters)

    R = initR()
    N = length(R[:,1]) #number of particles
    V = init_velocity(T,N) #velocity for t
    A = zeros(N,3,3) #acceleration matrix
    
    PE = Vector{Float64}(undef,iters+1) #initialize array of potential energies 
    # return fvector = update_A(R,T,V,dt) #initialize acceleration matrix and return initial PE

    PE[1], A = update_A(R,T,V,dt) #initialize acceleration matrix and return initial PE
    
    

    KE = Vector{Float64}(undef,iters+1) #initialize array of kinetic energies 
    KE[1], = calc_KE(V,N) #initial KE 

    T_list = Vector{Float64}(undef,iters+1) #initialize array of kinetic energies 
    T_list[1] = T #initial temperature


    for i in 1:iters
        curr_A = dropdims(sum(copy(A),dims=3),dims=3) #A for current t
        update_R!(R,V,curr_A,dt) #new R updated to t+dt
        PE[i+1],A = update_A(R,T,V,dt) #new A updated in-place to t+dt, and returns PE into array 
        Ac = copy(A[:,:,1])
        Ar = copy(A[:,:,2])
        update_V!(V,curr_A,Ac,Ar,dt) #new V updated to t+dt
        KE[i+1],T_list[i+1] = calc_KE(V,N) #add KE
        # T_list[i+1] = calc_T(KE[i+1],kb,N) #add T
    end
    return  PE, KE, T_list, A, V, R 
end

PE, KE, T_list, A, V, R= LD(2.,0.003,10000)