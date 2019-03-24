module KineticEnergy

using TypesBasis
using UtilityFunctions

export buildkineticenergy

@doc raw"""
    Use: buildkineticenergy(basisfunc::Basis)

Description: Generate the Ne*N-basis x Ne*N-basis kinetic energy matrix for our contracted basis set. The following expression captures this:

```
$ T_{n,m}$
```

Input:
Results:



"""
function buildkineticenergy(natoms::Int,basisfunc::Array{GaussOrbitals},basis::Basis)

    numbasisfunc = basis.nbasisfunc; #number of basis functions
    basissize = natoms*basis.nbasisfunc; #overlap size
    ke = zeros(Float64,basissize,basissize); #overlap matrix Ne*Gaussians x Ne*Gaussians

    flatbasisfunc = flattenbasisfunc(natoms,numbasisfunc,basisfunc);
  
    for n=1:basissize
        for m=1:basissize
            basisn = flatbasisfunc[n];
            basism = flatbasisfunc[m];
            
            #Gaussian factors
            an = basisn.alphas;
            am = basism.alphas;

            #Contraction coeffs.
            cn = basisn.coefs;
            cm = basism.coefs;

            xn,yn,zn = basisn.ax, basisn.ay, basisn.az;
            xm,ym,zm = basism.ax, basism.ay, basism.az;

            #Number of primitive Gaussians
            primsn = length(an);
            primsm = length(am);
            
            for np=1:primsn
                for mp=1:primsm
                    #Get overlap integral
                    overlap = getgaussoverlap(np,mp,basisn,basism);

                    #Contraction coeff prod.
                    cnm = cn[np]*cm[mp];
                    
                    #Gaussian product theorem
                    p = an[np]+am[mp];
                    px = (an[np]*xn + am[mp]*xm) / p;
                    py = (an[np]*yn + am[mp]*ym) / p;
                    pz = (an[np]*zn + am[mp]*zm) / p;                    
                    
                    #Build K.E. based on notes
                    ke[n,m] += 3 * cnm * am[mp] * overlap;
                    ke[n,m] -= 2 * cnm * am[mp]^2 * ((px-xm)^2 + 1/(2*p)) * overlap;
                    ke[n,m] -= 2 * cnm * am[mp]^2 * ((py-ym)^2 + 1/(2*p)) * overlap;
                    ke[n,m] -= 2 * cnm * am[mp]^2 * ((pz-zm)^2 + 1/(2*p)) * overlap;
                end          
            end # nprims        
        end 
    end #basissize    
    return ke
end #buildkineticenergy

end #KineticEnergy
