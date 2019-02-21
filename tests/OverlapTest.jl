push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis
using BuildBasis
using BuildOverlap
using InputParser

"""
		run_buildoverlap_test()

This is a function to test the user types in BuildOverlap.jl. In asserts overlap integral calculationfor two hydrogen atoms is the same as the output from Prof. James John's matlab code.

""" function run_buildoverlap_test()
    @testset "Overlap" begin
        system,basisset = getcalcsetup("../examples/HatreeFock.in")
        basisfunc,numelec = buildbasisfunc(system,basisset);
        Scalc = buildelecoverlap(numelec,basisfunc,basisset)
        Stest=  [1.00000e0 0.65829e0 0.77339e0 0.60892e0;
                 0.65829e0 1.00000e0 0.60892e0 0.95331e0;
                 0.77339e0 0.60892e0 1.00000e0 0.65829e0;
                 0.60892e0 0.95331e0 0.65829e0 1.00000e0]
        #Loop over and check each entry
        tolerance = 1.0e-5;
        for i=1:4
            for j=1:4
                @test begin
                    sign(Scalc[i,j]) == sign(Stest[i,j]);
                    diff = Scalc[i,j]-Stest[i,j];
                    diff <= tolerance;
                end
	    end
        end
        
    end #testset
end #run_overlap_test()

	   
