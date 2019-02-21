
######## REWRITE ################

push!(LOAD_PATH,"../src")

using Test
using TypesBasis
using ReadBasis

"""
		run_readbasis_test()

This is a function to test reading modified EMSL NWchem basis file. 

Required File: ../examples/basissets/H.6-31G.mod

""" function run_readbasis_test()
	 @testset "Read Basis" begin
	     #H Gaussian 6-31
	     format="modnwchem";
	     filename="../examples/basissets/H.6-31G.mod";
	     basisref = [["     18.7311370              0.03349460       ",
	 	         "      2.8253937              0.23472695       ",
		         "      0.6401217              0.81375733       "],
		         ["      0.1612778              1.0000000        "]];
             nbasisfunc, basisfunc = readgaussparams(filename,format=format);

             @test size(basisref) == size(basisfunc)

             stringbasisref = filter(!isspace,join(vcat(basisref...)));
             stringbasisfunc = filter(!isspace,join(vcat(basisfunc...)));
             @test stringbasisref == stringbasisfunc;
             
	end #testset
end

"""
		run_readbasis_test()

This is a function to test the constructed gaussian basis datatype

Required File: ../examples/basissets/H.6-31G.mod


""" function run_retreadbasis_test()
    @testset "Return Read Basis" begin
        #H Gaussian 6-31
        basistype = "6-31"
        numprims = 4;
        format="modnwchem";
        filename="../examples/basissets/H.6-31G.mod"
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
        nbasisfunc = gaussbasis.nbasisfunc;
        basisfunc = gaussbasis.basisfunc;

        @test nbasisfunc == 2;
        
        @testset "Basis-Func. 1" begin
            @test basisfunc[1].info == basistype;
            @test basisfunc[1].nprims == 3;
            @test basisfunc[1].coefs == refcoef[1:3];
            @test basisfunc[1].alphas == refexpo[1:3];
         end
   
        @testset "Basis-Func. 2" begin
            @test basisfunc[2].info == basistype;
            @test basisfunc[2].nprims == 1;
            @test basisfunc[2].coefs[1] == refcoef[4];
            @test basisfunc[2].alphas[1] == refexpo[4];   
        end
        
        
    end #testset
end #run_retreadbasis_test

    
            
