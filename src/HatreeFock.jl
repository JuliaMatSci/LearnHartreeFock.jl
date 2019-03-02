module HatreeFock

if  ! isempty(ENV["LEARNHATREEFOCK_PATH"])
    const SRC_PATH=ENV["LEARNHATREEFOCK_PATH"];
else
    println("! Environmental variable $LEARNHARTREEFOCK_PATH not set, defaulting to pwd() !")
    const SRC_PATH=pwd();
end

#using CodeInfo
include("ModuleList.jl")

"""
                        HatreeFock(filename=user_filename)

Description: This is the main function responsible for executing the Hatree-Fock
routines for calculating the total energy of a given system.

Program Outline
1.  Parse atomic system info:
             a.) parse(filename) returns a datatype with atomic number and coordinates.
2.  Caluclate nuclear-nuclear repulsion
3.  Setup basis functions:
            a.) returns the basis-set functions and number of electrons
            b.) build orbital overlap given to handle non-ortho. Gaussian basis set
4.   Calculate kinetic energy of each electron
5.   Build/calculate electron-nuclear attraction
6.   Build/calculate electron-electron repulsion (mean-field and exchange)
7.   Self-consistent field solution for total energy

""" function main(filename="HatreeFock.in")

    #printcodeinfo()
    
    atomicsystem, basis = getcalcsetup(filename);

    ZZ = calcnuclrepul(atomicsystem);

    basisfunc, numelec = buildbasisfunc(atomicsystem,basis);

    S = buildelecoverlap(numelec,basisfunc,basis);

    KE = buildkineticenergy(numelec,basisfunc,basis);

    Zq = buildnuclrattract(atomicsystem,numelec,
                             basisfunc,basis);

    qq = buildelecelecrepulsion(numelec,basisfunc,basis);

    size(KE) == size(Zq) ? nothing : throw(AssertionError("Kinetic energy matrix size does not equal nuclear attraction matrix!"))
    
    Ho = KE + Zq ;
    
    energyscf = runscf(numelec,S,Ho,qq);

    E = energyscf + ZZ;

    println("SCF Minimum energy: $E [Ha]")

    return 0
end #main

#Run Program
   main(ARGS[1])
#End Program

end
