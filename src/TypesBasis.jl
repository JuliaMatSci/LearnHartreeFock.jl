module TypesBasis

export Basis,BasisFunctions
export Gaussian,GaussOrbitals

using TypesParticles

abstract type Basis end

"""
Datatype for storing an array of basis functions.
""" struct BasisFunctions <: Basis
    nbasisfunc::Int
    basisfunc::Array{Basis,1}
    function BasisFunctions(nbasisfunc,basisfunc)
        if nbasisfunc != length(basisfunc)
            error("Length of basis function array does not equal the number of identified basis functions.")
        else
            new(nbasisfunc,basisfunc)
        end
    end
end #BasisFunction

        
""" 
Contains basis set info for Gaussian Functions

""" struct Gaussian <: Basis
       info::String
       nprims::Int
       coefs::Array{Float64,1}
       alphas::Array{Float64,1}
       function Gaussian(btype,nprims,coefs,alphas)
           if nprims != length(alphas)
               error("Number of alphas not equal to number of primitives.")
           elseif length(coefs) != length(alphas)
               error("Number of contractions not equal to number of alphas.")
           end
           new(btype,nprims,coefs,alphas)
       end
end #Gaussian



"""
Data type for Gaussian based orbitals centered on  ax, ay, az. No contractions or 
polarization functiosn just simple function form:

\$\\phi(x) = A e^{-\\alpha\\left(x-x_o\\right)^{2}}\$

""" struct GaussOrbitals
     ax::Float64
     ay::Float64
     az::Float64
     coefs::Array{Float64,1}
     alphas::Array{Float64,1}
     norms::Array{Float64,1}
     function GaussOrbitals(ax,ay,az,basis)
          norms = (2.00e0.*basis.alphas./pi).^(3/4);
          new(ax,ay,az,basis.coefs,basis.alphas,norms);
      end
     function GaussOrbitals(atom::Atom,basis)
          norms = (2.00e0.* basis.alphas ./pi).^(3/4);
          new(atom.x,atom.y,atom.z,basis.coefs,basis.alphas,norms);
      end
end #GausSOrbitals


end #module
