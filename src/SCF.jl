module SCF

#Julia Libraries/Modules
using Printf
using LinearAlgebra

#LearnHatreeFock.jl Modules
using TypesBasis,TypesParticles
using CoulombExchange

export runscf

@doc raw"""
function scf()

description: Self-consistent field soultion approach.

The Hatree-Fock operator for non-interacting electrons in an orthonormal basis set is give as:

$H_o C = \epsilon S C$ where $C,epsilon$ are the eigenvalues,eigenvectors and $S$ is the overlap matrix.

The first step is to find the eigenvalues and vectors of the overlap matrix. This will be used to transform from
atomic centered basis-set to a orthonormal molecular basis-set, for example, $H_{MO}=
"""
function runscf(numelec::Int,overlap::Array{Float64,2},
             ho::Array{Float64,2},elecelecrepul::Array{Float64,4};
             maxcycle=30::Int,tolerance=1.0e-6::Float64)

    convergflag = false;

    #Orthogonal overlap matrix, i.e., molecular orbitals
    oevals,oevecs = eigen(overlap);
    oevalsmat = Diagonal(oevals);

    #Try to invert eigenvalues matrix assuming non-singular, if singular
    # promote datatype ComplexF64
    invoevalsmat = returnevalsinvhalf(oevalsmat);
    orthoverlap = oevecs*invoevalsmat*adjoint(oevecs);
    orthoverlapstar = adjoint(orthoverlap);

    #Fock matrix construction
    fock = copy(ho);
    fockprime = orthoverlapstar*fock*orthoverlap;

    #HF Roothaan eq.
    #Find molecular orbital single electron solutions
    evals,evecsprime = eigen(fockprime);

    #Go back to atomic orbitals
    evecs = orthoverlap*evecsprime;
    
    #Sort lowest eigen values of for atomic orbitals
    sorteigen!(evals,evecs);

    
    #Build electron density for Hatree potential?
    density = builddensity(evecs,numelec)

    
    totalenergy = zeros(Float64,maxcycle+1);
    #SCF cycle
    for cycle=2:maxcycle+1

        #Get exchange and fock potential 
        exchangepot = getcoulexchange(density,elecelecrepul)
        fock = ho + exchangepot;

        #find solutions
        fockprime = orthoverlapstar*fock*orthoverlap;
        evals,evecsprime = eigen(fockprime);
        evecs = orthoverlap*evecsprime;
        
        sorteigen!(evals,evecs);
        
        density = builddensity(evecs,numelec);

        totalenergy[cycle] = getfockenergy(ho,fock,density);

        #Check conditions
        delta = abs(totalenergy[cycle]-totalenergy[cycle-1]);
        printscf(cycle-1,delta,totalenergy[cycle])
        if delta < tolerance
            return totalenergy[cycle]
        end
        
    end #cycle=2:maxcycle+1
    
    return totalenergy[end];
        
end #runscf


@doc raw"""
function returnevalsinvhalf(evals)

Try to invert eigenvalues matrix assuming non-singular, if singular promote 
datatype ComplexF64

"""
function returnevalsinvhalf(evals)
    try
        invevalsmat = evals^-0.5e0;
    catch
        invevalsmat = (convert(Array{ComplexF64},evals))^(-0.5e0);
    end
end #returnevalsinvhalf


@doc raw"""
function builddensity(evecs::Array{Number,2},n::)

Restricted HF density

"""
function builddensity(evecs::Matrix,numelec::Int)
    #Restricted HF density
    nhalf = Int(numelec/2);
    esize = size(evecs)[1];
    
    if typeof(evecs) == Array{Complex,2}
        density = zeros(Complex,esize,esize)
    else
        density = zeros(esize,esize)
    end
    for n=1:esize
        for m=1:esize
            for i=1:nhalf
                density[n,m] += evecs[n,i]*evecs[m,i];
            end #i=1:nhalf
        end
    end #n=1:esize
    return density
    
end #builddensity

@doc raw"""
function sortevals!(evals::Array,numelec::Int)

description: sort the eignevalues and vectors.

NOTES: I would like to rework this so that its more 
    succinct.
"""
function sorteigen!(evals::Vector{T},evecs::Matrix{T}) where {T<:Real}
    ncol = size(evecs)[2];

    #Make a shallow copy and force local scope
    #CHECK: I'm not sure what the most Julianic way is to do this.
    local sortedevals = copy(evals);
    local sortedevecs = copy(evecs);
    
    #Get sorted indexes for eigenvals
    sortedindex = sortperm(evals);
    evals[:] = sortedevals[sortedindex];

    #Get sorted eigenvectors
    for i=1:ncol
        sortedevecs[:,i] = evecs[:,sortedindex[i]];
    end

    evecs[:,:] = sortedevecs[:,:];

end #sortevals

function sorteigen!(evals::Vector{T},evecs::Matrix{T}) where {T<:Complex}
    ncol = size(evecs)[2];

    #Make a shallow copy and force local scope
    #CHECK: I'm not sure what the most Julianic way is to do this.
    local sortedevals = copy(evals);
    local sortedevecs = copy(evecs);

    magevals = adjoint(evals)*evals;
    println(magevals);
    #Get sorted indexes for eigenvals
    sortedindex = sortperm(real(magevals));
    evals[:] = sortedevals[sortedindex];

    #Get sorted eigenvectors
    for i=1:ncol
        sortedevecs[:,i] = evecs[:,sortedindex[i]];
    end

    evecs[:,:] = sortedevecs[:,:];

end #sortevals

@doc raw"""
sorteigen()

TODO
"""
function sorteigen(evals::Vector{T},evecs::Matrix{T}) where {T<:Real}
    p = sortperm(evals);
    return evals[p], evecs[p,:]
end

@doc raw"""
function getfockenergy()

"""
function getfockenergy(ho::Array,fock::Array,density::Array)
    numbasis = size(density)[1];
    energy = 0.00e0;
    for n=1:numbasis
        for m=1:numbasis
            energy += density[n,m] * (ho[n,m]+fock[n,m]);
        end
    end#n=1:numbasis
    return energy    
end #getfockenergy

@doc raw"""
function printscf()

"""
function printscf(iteration::Int,delta::Float64,energy::Float64)
    if iteration < 2
        @printf("                   HF SCF                         \n");
        @printf("--------------------------------------------------\n");
        @printf("Iteration         ΔE [Ha]          HF-Energy [Ha] \n");
        @printf("--------------------------------------------------\n");
        @printf("%i                %0.6f           %0.6e   \n",
            iteration,delta,energy);
    else
        @printf("%i                %0.6f           %0.6e   \n",
            iteration,delta,energy);
    end
end #printscf

end #SCF
