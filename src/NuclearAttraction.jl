module NuclearAttraction

#Julia external modules
using GSL

using TypesBasis,TypesParticles
using UtilityFunctions

export buildnuclrattract

@doc raw"""
function buildnuclrattract(system::)

This function applies the Coulomb operator between nuclei and electronic volume, i.e., the volume
integral of the electron density.


```
$V = int d\tau \frac{-Z_{i}}{r_{i,n}\phi^2_{n}(r_{j})$
```

This function will only support atom centered basis functions and thus the potential is:

```
V_{i,j} = \sum_N - Z_n \int\frac{\phi_i \phi_j}{r_{i,j}} d\tau

```

The overlap $phi_i \phi_j$ can be written as $\phi_p = N_p e^{-pr_{p}^2} and the integral becomes
$\int\frac{\phi_p}{r_{i,j}}$

The approach to solving this for spatial integral of gaussian functions is to use Boys functions:

```
#\int \frac{\phi_p}{r_{i,j}} = \frac{2\pi C_p}{p} F_o(pR^2_{NP})
```

Input:
Results:

module requirements:

"""
function buildnuclrattract(system::SystemOfAtoms,basisfunc::Array{GaussOrbitals},basis::Basis)
    natoms = system.natoms;
    numbasisfunc = basis.nbasisfunc; #number of basis functions
    basissize = natoms*numbasisfunc;
    nuclrattract = zeros(Float64,basissize,basissize);
    
    flatbasisfunc = flattenbasisfunc(natoms,numbasisfunc,basisfunc);

    for i=1:natoms
        atom = system.atoms[i];
        for n=1:basissize
            for m=1:basissize
                
                #Basis functions with primitives
                basisn = flatbasisfunc[n];
                basism = flatbasisfunc[m];

                #Contraction coeffs.
                cn = basisn.coefs;
                cm = basism.coefs;

                #Number of primitive Gaussians
                primsn = length(cn);
                primsm = length(cm);

                for np=1:primsn
                    for mp=1:primsm
                        cnm = cn[np]*cm[mp];
                        nuclrattract[n,m] += cnm * coulomb(np,mp,basisn,basism,atom);
                    end
                end #primsn

            end
        end #basissize
    end #system.natoms
    return nuclrattract
    
end #buildnuclrattract

@doc raw"""
function coulomb(np::Int,mp::Int,basisn::Gaussian
coulomb potential

"""
function coulomb(np::Int,mp::Int,basisn::GaussOrbitals,
                 basism::GaussOrbitals,atom::Atom)

    xn,yn,zn = basisn.ax,basisn.ay,basisn.az;
    xm,ym,zm = basism.ax,basism.ay,basism.az;

    alphanp = basisn.alphas[np];
    alphamp = basism.alphas[mp];
    
    normcont = basisn.norms[np] * basism.norms[mp];
    
    p = alphanp + alphamp;
    px = (alphanp*xn + alphamp*xm)/p;
    py = (alphanp*yn + alphamp*ym)/p;
    pz = (alphanp*zn + alphamp*zm)/p;
    orbatomdist = (atom.x-px)^2 + (atom.y-py)^2 + (atom.z-pz)^2;

    spatialx = gaussprod1D(xn,alphanp,xm,alphamp);
    spatialy = gaussprod1D(yn,alphanp,ym,alphamp);
    spatialz = gaussprod1D(zn,alphanp,zm,alphamp);

    spatial = spatialx*spatialy*spatialz;
    prefact = (-1*atom.Z*2*pi) / p;

    ##WARNING! call is for spherical Gaussians only, x=0.
    potential = prefact * normcont * spatial * boysfunction(p*orbatomdist,n=0);

    return potential
    
end #coulomb

end #NuclearAtraction
