module BuildBasis

export buildbasisfunc

using TypesParticles, TypesBasis


@doc raw"""
buildbasisfunction(atomicsystem::SystemOfAtoms,basis::Basis)

Returns atom-centered basis functions and number of electrons. 

Description: This function will take the system information and create the basis functions
that are atomically centered. Each basis function will contain the parameters for the Gaussian
functions that when linearly combined make up the basis function. Note that the actual basis 
function is not yet constructed here, i.e., we are not taking the summation of Gaussians. 

Depends on datatypes in TypesParticles.jl and TypesBasis.jl

""" function buildbasisfunc(atomicsystem::SystemOfAtoms,basis::Basis)
    natoms = atomicsystem.natoms;
    nelectrons = 0;

    #Closed-shell system only
    for i=1:natoms
        nelectrons += atomicsystem.atoms[i].Z
    end
    
    nbasis = basis.nbasisfunc;
    basisfunc = Array{GaussOrbitals,2}(undef,natoms,nbasis);

    #Assign each atom-centered basis functions
    #ie = 1; #total electron counter
    for a=1:natoms
        #For each gaussian store 5 values: c,alpha,xo,yo,zo
        #For Z=1 n-basis per atom is 1 so basisfunc length equal num atoms. 
        for j=1:nbasis
            #6-31 type basis set
            if basis.basisfunc[j].info == "6-31"
                basisfunc[a,j] =  GaussOrbitals(atomicsystem.atoms[a],basis.basisfunc[j])
            end
        end # nbasis
    end # atomicsystem.atoms
    
    return basisfunc,nelectrons
end #buildbasisfunction

end #BuildBasis
