module BuildBasis

export buildbasisfunc

using TypesParticles, TypesBasis


"""
                      build_basis_function(atomicsystem::SystemOfAtoms,basis::Basis)

returns atom-centered basis functions and number of electrons. 

Description: this function will take the system information and create the basis functions
that are atomically centered. Each basis function will contain the parameters for the Gaussian
functions that when linearly combined make up the basis function. Note that the actual basis 
function is not yet constructed here, i.e., we are not taking the summation of Gaussians. 

""" function buildbasisfunc(atomicsystem::SystemOfAtoms,basis::Basis)
    #basis is a BasisFunctions datatype
    
    natoms = atomicsystem.natoms;
    nelectrons = 0;
    for i=1:natoms
        nelectrons += atomicsystem.atoms[i].Z
    end
    
    nbasis = basis.nbasisfunc;
    basisfunc = Array{GaussOrbitals,2}(undef,nelectrons,nbasis);

    #Assign each atom-centered basis functions
    for i=1:nelectrons
        #Get correct expon.
        #For each gaussian store 5 values: c,alpha,xo,yo,zo

        #For Z=1 n-basis per atom is 1 so basisfunc length equal num atoms. 
        for j=1:nbasis
            #6-31 type basis set
            if basis.basisfunc[j].info == "6-31"
                basisfunc[i,j] =  GaussOrbitals(atomicsystem.atoms[i],basis.basisfunc[j])
            end
        end # nbasis
         
    end #nelectrons
    
    return basisfunc,nelectrons
end #build_basis_function

    
#"""
#Multiple dispatch version for Array of Basis type

#""" function buildbasisfunc(atomicsystem::SystemOfAtoms,basis::Array{T,1}) where T <: Basis
#    natoms = atomicsystem.natoms;
#    nelectrons = 0;
#    nbasis = 0;
#    nbasisorb = length(basis)    
#    #Assign each atom-centered basis functions
#    basisfunc = Array{GaussOrbitals,2}(undef,natoms,nbasisorb);
#     for i=1:natoms
#        Z=atomicsystem.atoms[i].Z
#        nelectrons += Z
#        #Get correct expon.
#        #For each gaussian store 5 values: c,alpha,xo,yo,zo

#        #For Z=1 n-basis per atom is 1 so basisfunc length equal num atoms. 
#        for n=1:nbasisorb
#            if (basis[n].info == "6-31") && ( Z == 1)
#                basisfunc[i,n] =  GaussOrbitals(atomicsystem.atoms[i],basis[n])
#            end # "6-31"           
#        end # nbasisorb
         
#    end
#    return basisfunc,nelectrons
#end #build_basis_function


end #module
