module NuclearRepulsion

using TypesParticles

export calcnuclrepul

"""
    Use: calcnuclrepul(system:SystemOfAtoms)

Input:
Results:

Description: returns the total energy from the nuclear-nuclear repulsion. This is a direct
sum calculation, as opposed to Ewald summation or particle mesh solvers.

TODO: return energy in units of eV or Ha

""" function calcnuclrepul(system::SystemOfAtoms)
    energy = 0.00e0;
    if system.natoms == 1
        return energy
    else
        for i=1:system.natoms
            atomi = system.atoms[i]
            for j=i+1:system.natoms;
                atomj = system.atoms[j];
                mag_rij = sqrt((atomj.x-atomi.x)^2 +
                              (atomj.y-atomi.y)^2 +
                              (atomj.z-atomi.z)^2);
                energy = energy + atomi.Z*atomj.Z/mag_rij;
            end
        end #system.natoms
        return energy
    end
end #get_nucl_nucl


end #NuclearRepulsion
