push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis
using BuildBasis
using InputParser
using ElecElecRepulsion

"""
		run_buildelecelecrepul_test()



""" function run_buildelecelecrepul_test()
    @testset "Electron-Electron Repulsion" begin
       system,basisset = getcalcsetup("../examples/HatreeFock.in")
       basisfunc,numelec = buildbasisfunc(system,basisset);
       elecelecrepcalc = buildelecelecrepulsion(numelec,basisfunc,basisset);
       @test 1 == 1
        # tolerance = 1.0e-4;
        # for i=1:4
        #     for j=1:4
        #         @test sign(nuclrattractcalc[i,j]) == sign(nuclrattracttest[i,j]);
        #         diff = abs(nuclrattractcalc[i,j]-nuclrattracttest[i,j]);
        #         @test diff <= tolerance;
	#     end
        # end #i

    end #testset
    println("!!!!! TEST NOT IMPLEMENTED YET !!!!!")
end #run_elecelecrepul_test()

	   
