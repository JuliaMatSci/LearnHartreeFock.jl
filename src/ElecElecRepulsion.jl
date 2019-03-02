module ElecElecRepulsion

#LearnHatreeFock.jl Modules
using TypesBasis,TypesParticles
using UtilityFunctions

export buildelecelecrepulsion

@doc raw"""
function buildelecelecrepulsion()

description: Implement electron-electron repulsion within Hatree-Fock framework. The electron-electron repuslion is done within a mean-field approach. That is given the electron density what is the average repulsion energy an electron feels due to the charge density located at a different point. Recall that charge density is given by the overlap of two orbitals. Our operator is thus:

``` 
$ \int e^2 frac{\psi_{m}}{\psi_{n}} dr $
```

We need to have all possible combinations of orbitals for the electron-electron operator, N^4.

Gaussian orbitals centered at different locations. Break-down ad two set of orbitals using gaussian product theorem, so we get two gaussian at new centers called p and p-prime. The integral is then given by:

```
$g(l,m,n,o) = A_p\,A_{p^{`}} frac{2\pi^{\frac{5}{2}}}{p\,p^{`}\,\sqrt{p+p^{`}}} F_o\left(\alpha_{p,p^{`}}\,R_{p,p^{`}}^2 \right)$        
```

The steps for solving the integral are:
1. Loop over all basis functions, i.e., l,m,n,o
2. Loop over all primitive Gaussian functions of l,m,n,o.
3. Apply Gaussian product theorem to p=lm and p`=no

input:
output:


"""
function buildelecelecrepulsion(numelec::Int,basisfunc::Array{GaussOrbitals},
                                basis::Basis)

    numbasisfunc = basis.nbasisfunc; #number of basis functions
    basissize = numelec*numbasisfunc;

    # Ne*Nbasis x Ne*Nbasis x Ne*Nbasis x Ne*Nbasis 
    elecelecrepul = zeros(Float64,basissize,basissize,
                          basissize,basissize);
    
    flatbasisfunc = flattenbasisfunc(numelec,numbasisfunc,basisfunc);

    #NOTE: To try and comp. cost make min. call/ref.
    #Build electron-electron Coulomb intergral
    for l=1:basissize
        basisl = flatbasisfunc[l];
        primsl = length(basisl.coefs);
        for lp=1:primsl
            for m=1:basissize
                basism = flatbasisfunc[m];
                primsm = length(basism.coefs);
                for mp=1:primsm
                    #Build new Gaussian center from l,m Gaussian orbitals.
                    #TODO: Turn into function
                    gausslm = getnewcenter(basisl,lp,basism,mp);

                    #Now build new Gaussian center for n,o
                    for n=1:basissize
                        basisn = flatbasisfunc[n];
                        primsn = length(basisn.coefs);
                        for np=1:primsn
                            for o=1:basissize
                                basiso = flatbasisfunc[o];
                                primso = length(basiso.coefs);
                                for op=1:primso

                                    gaussno = getnewcenter(basisn,np,basiso,op);
                                    alpha = (gausslm.p*gaussno.p)/(gausslm.p+gaussno.p);
                                    dist = (gausslm.px-gaussno.px)^2 +
                                           (gausslm.py-gaussno.py)^2 +
                                           (gausslm.pz-gaussno.pz)^2;
                                    #NOTE: This calls a spherical implemented version
                                    boysfunc = boysfunction(alpha*dist,n=0);
                                    prefact = (2.00e0*pi^2.50e0)/(gausslm.p*gaussno.p*sqrt(gausslm.p+gaussno.p));
                                    elecelecrepul[l,m,n,o] += prefact*gausslm.ampli*gaussno.ampli*boysfunc;
                                end
                            end # o=1:basissize
                        end
                    end # n=1:basissize
                end
            end # m=1:basissize
        end
    end # l=1:basissize
    
    return elecelecrepul
end #buildelecelecrepulsion


@doc raw"""
function getnewcenter();

"""
function getnewcenter(basisn::GaussOrbitals,np::Int,
                      basism::GaussOrbitals,mp::Int)
    alpha_np = basisn.alphas[np]::Float64;
    alpha_mp = basism.alphas[mp]::Float64;
    coef_np = basisn.coefs[np]::Float64;
    coef_mp = basism.coefs[mp]::Float64;
    norm_np = basisn.norms[np]::Float64;
    norm_mp = basism.norms[mp]::Float64;
    p = alpha_np+alpha_mp;
    px = (basisn.ax*alpha_np + basism.ax*alpha_mp)/p;
    py = (basisn.ay*alpha_np + basism.ay*alpha_mp)/p;
    pz = (basisn.az*alpha_np + basism.az*alpha_mp)/p;
    spatialx_nm =gaussprod1D(basisn.ax,alpha_np,basism.ax,alpha_mp);
    spatialy_nm =gaussprod1D(basisn.ay,alpha_np,basism.ay,alpha_mp);
    spatialz_nm =gaussprod1D(basisn.az,alpha_np,basism.az,alpha_mp);
    spatial_nm = spatialx_nm*spatialy_nm*spatialz_nm;
    ampli_nm = norm_np * norm_mp * coef_np * coef_mp * spatial_nm;
    return NewGauss(p,px,py,pz,ampli_nm)

end #getnewcenter

"""
TODO: should I move this to BasisTypes.jl?

""" struct NewGauss
     p::Float64
     px::Float64
     py::Float64
     pz::Float64
     ampli::Float64
     function NewGauss(p,px,py,pz,ampli)
          new(p,px,py,pz,ampli);
      end
end #NewGauss

end #ElecElecRepulsion
