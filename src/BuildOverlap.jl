module BuildOverlap

export buildelecoverlap

using TypesParticles, TypesBasis
using  UtilityFunctions

@doc raw"""
                      buildelecoverlap(basis::Basis)

returns the building of electron overlap integral matrix of expectation values:

     S_{m,n} = langle m | n rangle = int d tau phi_m phi_n 

Description: This function calculates the electron-electron overlap matrix which essential captures
the inner product projection of basis function m on n, this is neccessary to track basis-functions that
are not orthogonal to one another.

NOTES: multiple dispatch only for GaussOrbitals

"""
function buildelecoverlap(numelec::Int,basisfunc::Array{GaussOrbitals},basis::Basis)
    numbasisfunc = basis.nbasisfunc; #number of basis functions
    osize = numelec*basis.nbasisfunc; #overlap size
    overlap = zeros(osize,osize); #overlap matrix Ne*Gaussians x Ne*Gaussians

    #flatbasisfunc=vcat(basisfunc...); #Need to improve
    flatbasisfunc = flattenbasisfunc(numelec,numbasisfunc,basisfunc);
  
    for n=1:osize
        for m=1:osize
            nbasis = flatbasisfunc[n];
            mbasis = flatbasisfunc[m];
            nprims = length(nbasis.alphas);
            mprims = length(mbasis.alphas);    
            for np=1:nprims
                for mp=1:mprims
                    coefs = nbasis.coefs[np]*mbasis.coefs[mp];
                    gaussoverlap = getgaussoverlap(np,mp,nbasis,mbasis);
                    overlap[n,m] += coefs*gaussoverlap;
                end
            end # mprims        
        end 
    end #osize
    
    return overlap
end #buildelecoverlap

end #module
