push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis
using BuildBasis
using BuildOverlap
using InputParser
using KineticEnergy

"""
		run_buildkineticenergy_test()

This is a function to test the user types in BuildOverlap.jl. In asserts overlap integral calculationfor two hydrogen atoms is the same as the output from Prof. James John's matlab code.

""" function run_buildkineticenergy_test()
    @testset "Kinetic Energy" begin
        system,basisset = getcalcsetup("../examples/HatreeFock.in")
        basisfunc,numelec = buildbasisfunc(system,basisset);

        Scalc = buildelecoverlap(numelec,basisfunc,basisset);
        kecalc = buildkineticenergy(numelec,basisfunc,basisset);
        ketest = [1.39568e0 0.25974e0 0.75367e0 0.22767e0;
                  0.25974e0 0.24192e0 0.22767e0 0.22327e0;
                  0.75367e0 0.22767e0 1.39568e0 0.25974e0;
                  0.22767e0 0.22327e0 0.25974e0 0.24192e0];

        tolerance = 1.0e-5;
        for i=1:4
            for j=1:4
                @test sign(kecalc[i,j]) == sign(ketest[i,j]);
                diff = abs(kecalc[i,j]-ketest[i,j]);
                @test diff <= tolerance;
	    end
        end #i
        
    end #testset
end #run_kineticenergy_test()

	   
