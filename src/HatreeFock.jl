module HatreeFock

using CodeInfo

"""
                        HatreeFock(filename=user_filename)

Description: This is the main function responsible for executing the Hatree-Fock
routines for calculating the total energy of a given system.

Program Outline
1.  Parse atomic system informatics:
             a.) parse(filename) returns a datatype with atomic number and coordinates.
2.  Caluclate nuclear-nuclear repulsion
3.  Build basis overlap:
            a.) returns the basis-set functions and number of electrons
4.   Calculate kinetic energy of each electron
5.   Build/calculate electron-nuclear attraction
6.   Build/calculate electron-electron repulsion (mean-field and exchange)
7.   Self-consistent field solution for total energy

""" function main(filename="HatreeFock.in"))

    printcodeinfo()
    
    atomicsystem, basis = getcalcsetup(filename);

    ZZ = calcnuclrepul(atomicsystem);

    basisfunc, numelec = build_basis_func(atomicsystem,basis);

    S = buildelecoverlap(basisfunc,basis);

    KE = build_kinetic_energy(basisfunc);

    Zq = build_nuclear_attract(atomicsystem,basisfunc);

    qq = build_elec_elec_repulsion(basisfunc);

    H_o = KE + ZZ ;
    
    E_scf = min_energy_via_scf(H_o,qq,S,numelec);

    E = E_scf + ZZ;

    println("SCF Minimum energy: $E [Ha]")

    return 0
end #main

#Run Program
   main(ARGS[1])
#End Program

end
