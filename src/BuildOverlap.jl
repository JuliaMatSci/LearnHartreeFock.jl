module BuildOverlap

export buildelecoverlap

using TypesParticles, TypesBasis

"""
                      buildelecoverlap(basis::Basis)

returns the building of electron overlap integral matrix of expectation values:

     S_{m,n} = langle m | n rangle = int d tau phi_m phi_n 

Description: This function calculates the electron-electron overlap matrix which essential captures
the inner product projection of basis function m on n, this is neccessary to track basis-functions that
are not orthogonal to one another.

NOTES: multiple dispatch only for GaussOrbitals

""" function buildelecoverlap(numelec::Int,basisfunc::Array{GaussOrbitals},basis::Basis)
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

"""

TODO write description

""" function getgaussoverlap(np::Int,mp::Int,
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
        
"""
TODO write description

""" function gaussprod1D(npx::Float64,np_alpha::Float64,
                         mpx::Float64,mp_alpha::Float64)
    p = np_alpha + mp_alpha;
    q = (np_alpha * mp_alpha) / p;
    #ratio = (np_alpha*npx + mp_alpha*mpx)/ p;
    diff = npx-mpx;
    kernel = exp(-q*diff^2);
    return kernel
end #gaussprod1D


"""
TODO write description

""" function flattenbasisfunc(numelec::Int,numbasisfunc::Int,
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


end #module
