module UtilityFunctions

export gaussprod1D,getgaussoverlap
export flattenbasisfunc
export boysfunction

using TypesParticles, TypesBasis

@doc raw"""
TODO write description

"""
function gaussprod1D(npx::Float64,np_alpha::Float64,
                         mpx::Float64,mp_alpha::Float64)
    p = np_alpha + mp_alpha;
    q = (np_alpha * mp_alpha) / p;
    #ratio = (np_alpha*npx + mp_alpha*mpx)/ p;
    diff = npx-mpx;
    kernel = exp(-q*diff^2);
    return kernel
end #gaussprod1D

@doc raw"""

TODO write description

"""
function getgaussoverlap(np::Int,mp::Int,
                             nbasis::GaussOrbitals,mbasis::GaussOrbitals)
    np_alpha = nbasis.alphas[np];
    mp_alpha = mbasis.alphas[mp];
    spatialx =  gaussprod1D(nbasis.ax,np_alpha,mbasis.ax,mp_alpha);
    spatialy =  gaussprod1D(nbasis.ay,np_alpha,mbasis.ay,mp_alpha);
    spatialz =  gaussprod1D(nbasis.az,np_alpha,mbasis.az,mp_alpha);
    sprod = spatialx*spatialy*spatialz;
    term1 = sqrt(pi/(np_alpha+mp_alpha))^3; #relabel based on equations
    term2 = nbasis.norms[np]*mbasis.norms[mp]; #rekabek based on equations
    primoverlap = sprod*term1*term2;
    return primoverlap
end #getgaussoverlap
        

@doc """
flattenbasisfunction(number_electrons::Int, number_basis_functions::Int,
                     basis_functions::Array{GaussOrbitals})

description: This utility function takes a 2-dimensional array with elements having data-type
{GaussOrbitals} and flattens it to a 1-dimensional of data-type {Gaussorbitals} that has as
length equal to number of electrons times number of basis functions.

input:

output:

module requirements: TypesParticles.jl, TypesBasis.jl

"""
function flattenbasisfunc(numelec::Int,numbasisfunc::Int,
                              basisfunc::Array{GaussOrbitals})
    overlapsize = numelec*numbasisfunc;
    flatbasisfunc = Array{GaussOrbitals,1}(undef,overlapsize);
    ei = 1;
    for e=1:numelec
        tmpbf = basisfunc[e,:]
        for i=1:numbasisfunc
            flatbasisfunc[ei] = tmpbf[i];
            ei += 1;
        end
    end
    return flatbasisfunc
end #flattenbasisfunc

@doc raw"""
function boysfunction(x::Float64;n::Int=0)

description: This is an implementation for integration of coulomb basis intergral using
boys functions for Gaussian,  that are atomic-centered, basis functions. Requires use of gamma and incomplete gamma functions obtained through GSL library..

input: 
output:

module requirements: This function needs access to the gamma and incomplete gamma functions. For this function call I use the GSL (GNU Scientific Library) interface which has sf_gamma and complimentary gamma function sf_gamma_inc_P

"""
function boysfunction(x::Float64;n::Int=0)
    if x == 0.00e0
        out = 1.00e0/(2.00e0*n+1.00e0);
    else
        nhalf = n + 0.5e0;
        gamma = GSL.sf_gamma(nhalf);
        igamma = GSL.sf_gamma_inc_P(nhalf,x);
        out = gamma*igamma/(2*x^(n+0.5e0));
    end
    return out
    
end #boysfunction


end #module
