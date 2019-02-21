push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis
using BuildBasis
using BuildOverlap
using InputParser

"""
		run_buildoverlap_test()

This is a function to test the user types in BuildOverlap.jl. In asserts overlap integral calculation is
the same as James John's matlab code

""" function run_buildoverlap_test()
#    @testset "Overlap" begin
#        @test begin
            system,basisset = getcalcsetup("./HatreeFock.in")
            basisfunc,numelec = buildbasisfunc(system,basisset);
            S = buildelecoverlap(numelec,basisfunc,basisset)
            println(S);
            Stest=  [1.00000   0.65829   0.77339   0.60892;
                     0.65829   1.00000   0.60892   0.95331;
                     0.77339   0.60892   1.00000   0.65829;
                     0.60892   0.95331   0.65829   1.00000]
            Scalc = [1.0 0.658292 0.77339 0.608916;
                     0.658292 1.0 0.608916 0.953314;
                     0.77339 0.608916 1.0 0.658292;
                     0.608916 0.953314 0.658292 1.0]
    
#	end
#    end #testset
end #run_overlap_test()

	   
