push!(LOAD_PATH,"../src")

using Test
using TypesParticles
using TypesBasis

"""
		run_typesparticles_test()

This is a function to test the user types in TypesParticles.jl. In asserts basics
construction of the user data type for hydrogen.

""" function run_typesparticles_test()
	 @testset "Types" begin
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
	       	H1.x == x0;
	       	H1.y == y0;
	       	H1.z == z0;
	       	H1.Z == Z0;
	       	end
		
		@test begin
	       	atomlist = [H1,H2];
	       	system = SystemOfAtoms(2,atomlist,cell);
	       	system.natoms == 2;	
	       	typeof(system.atoms) == typeof(atomlist);
		system.atoms[1] == H1;
		system.atoms[2] == H2;
	       	system.cell == cell;
	       	end
		
	 end #testset
end #run_typeparticles_test


"""
		run_typesbasis_test()

This is a function to test the user types in TypesBasis.jl. In asserts basic
construction of the user data type for H 6-31G (Gaussian) basis set.

""" function run_typesbasis_test()
	 @testset "Basis" begin
	 #Hydrogen 6-31G
	 basistype="6-31G";
	 nprimatives=4;
	 coefficients = [ 0.03349460;
	 	          0.23472695;
			  0.81375733;
			  1.00000000];
	 alphas = [18.7311370;
	 	   2.8253937;
		   0.6401217;
		   0.1612778];
         norms = (2.00e0.*alphas./pi).^(3/4);
	 @test begin
	      basisset = Gaussian(basistype,nprimatives,
		       coefficients,alphas);
		       
	      basisset.info == basistype;
	      basisset.nprims == nprimatives;
	      basisset.coefs == coefficients;
	      basisset.alphas == alphas;
	 end

	 @test begin
	       x0,y0,z0 = 0.00e0, 0.00e0, 0.00e0;
	       Z0 = 1;
	       H1 = Atom(x0,y0,z0,Z0);
	       
	       basisset = Gaussian(basistype,nprimatives,
		       coefficients,alphas);
	       gaussorb1 = GaussOrbitals(x0,y0,z0,basisset);
	       gaussorb2 = GaussOrbitals(H1,basisset);

	       gaussorb1.ax == x0;
	       gaussorb2.ax == x0;
	       gaussorb1.ay == y0;
	       gaussorb2.ay == y0;
	       gaussorb1.az == z0;
	       gaussorb2.az == z0;
	       
	       gaussorb1.coefs == coefficients;
	       gaussorb2.coefs == coefficients;

	       gaussorb1.alphas == alphas;
	       gaussorb2.alphas == alphas;

	       gaussorb1.norms == norms;
	       gaussorb2.norms == norms;
	   end
	   end #testset
end #run_typesbasis_test()

	   
