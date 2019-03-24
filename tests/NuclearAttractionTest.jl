push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis
using BuildBasis
using InputParser
using NuclearAttraction

"""
		run_buildnuclrattract_test()



""" function run_buildnuclrattract_test()
    @testset "Nuclear Attraction" begin
        system,basisset = getcalcsetup("../examples/HatreeFock.in")
        basisfunc,numelec = buildbasisfunc(system,basisset);
        nuclrattractcalc = buildnuclrattract(system,basisfunc,basisset);
        nuclrattracttest = [-2.7494e0 -1.3415e0 -2.1471e0 -1.2642e0;
                  -1.3415e0 -1.2431e0 -1.2642e0 -1.2027e0;
                  -2.1471e0 -1.2642e0 -2.7494e0 -1.3415e0;
                  -1.2642e0 -1.2027e0 -1.3415e0 -1.2431e0];
        
        tolerance = 1.0e-4;
        for i=1:4
            for j=1:4
                @test sign(nuclrattractcalc[i,j]) == sign(nuclrattracttest[i,j]);
                diff = abs(nuclrattractcalc[i,j]-nuclrattracttest[i,j]);
                @test diff <= tolerance;
	    end
        end #i
        
    end #testset
end #run_nuclrattract_test()

	   
