module TypesParticles

export Particle, Atom
export SystemOfAtoms

abstract type Particle end

"""
ab-initio Information about an atom.

Use: Atom(x::Float64,y::Float64,z::Float64,Z::Integer)
""" struct Atom <: Particle
           x::Float64
           y::Float64
           z::Float64
           Z::Integer
           function Atom(x::Real=0.00e0,
                         y::Real=0.00e0,
                         z::Real=0.00e0,
                         Z::Real=1)
                   x0,y0,z0 = map(Float64,(x,y,z));
                   Z0 = Int(Z);
                   if Z0 < 1
                       error("Atomic number < 1 is not allowed.")
                   end
                   new(x0,y0,z0,Z0);
           end
end #Atom         

"""
Data type for storing list of information about an atomic system.

Use: SystemOfAtoms(natoms::Integer,atoms::Array{Atom,1},cell::Array{Float64,2}(3x3))
""" struct SystemOfAtoms
    natoms::Integer;
    atoms::Array{Atom,1};
    cell::Array{Float64,2};
    function SystemOfAtoms(natoms=1,
                    atoms=[Atom()],
                    cell=zeros(3,3))
            if natoms < 0
                error("Simulation cell cannot have negative atoms.")
            end
            if size(cell) != (3,3)
                error("Simulation cell is not a 3x3 array.")
            end
            if typeof(atoms) != Array{Atom,1}
                error("Atomic entries are not of the typeof() == Atom.")
            end
            new(natoms,atoms,cell)
    end
end #SystemOfAtoms

end#module


