push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using NuclearRepulsion

"""
	run_nuclearrepulsion()	

""" function run_nuclearrepulsion_test()
	 @testset "NuclearRepulsion" begin
	 #H2 structure
	 x0,y0,z0 = 0.00e0,0.00e0,0.00e0;
	 x1,y1,z1 = 0.77e0,0.00e0,0.00e0;
	 Z0,Z1 = 1,1;
	 H1 = Atom(x0,y0,z0,Z0);
	 H2 = Atom(x1,y1,z1,Z1);

	 cell = [1.00e0 0.00e0 0.00e0;
	      	 0.00e0 1.00e0 0.00e0;
		 0.00e0 0.00e0 1.00e0];
            
             @test begin
	         atomlist = [H1,H2];
	         system = SystemOfAtoms(2,atomlist,cell);
	         calcnuclrepul(system) == 1.2987012987012987e0; #Atomic Units
                 1 == 1
             end
        end #testset
end #run_nuclearrepulsion_test
