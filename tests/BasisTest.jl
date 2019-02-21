
######## REWRITE ################

push!(LOAD_PATH,"../src")

using Test
using TypesBasis
using ReadBasis

"""
		run_readbasis_test()

This is a function to test reading modified EMSL NWchem basis file. 

File: ./basisfiles/H.6-31G.mod

""" function run_readbasis_test()
	 @testset "Read Basis" begin
	     #H Gaussian 6-31
	     format="modnwchem";
	     filename="./basisfiles/H.6-31G.mod";
	     basisref = ["     18.7311370              0.03349460       ",
	 	         "      2.8253937              0.23472695       ",
		         "      0.6401217              0.81375733       ",
		         "      0.1612778              1.0000000        "];
             basisread = readgaussparams(filename,format=format);
	     @test isequal(basisref,basisread)	  
	end #testset
end

"""

""" function run_retreadbasis_test()
    @testset "Return Read Basis" begin
        #H Gaussian 6-31
        basistype = "6-31"
        numprims = 4;
        format="modnwchem";
        filename="./basisfiles/H.6-31G.mod"
        refexpo = [18.7311370,
                  2.8253937,
                  0.6401217,
                  0.1612778];
        refcoef = [0.03349460,
	 	   0.23472695,
		   0.81375733,
		   1.0000000];
        
        gaussbasis = retreadbasis(filename,format=format,
                                  btype=basistype)
        @test begin
            gaussbasis.info == basistype;
            gaussbasis.nprims == 4;
            gaussbasis.coefs == refcoef;
            gaussbasis.alphas == refexpo;
        end
    end #testset
end #run_retreadbasis_test

    
            
