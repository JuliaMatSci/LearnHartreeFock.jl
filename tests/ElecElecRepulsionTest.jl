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
       elecelecreptest = gettestdata();
       tolerance = 1.0e-4;
       for i=1:4
           for j=1:4
               for k=1:4
                   for l=1:4
                       @test sign(elecelecrepcalc[i,j,k,l]) == sign(elecelecreptest[i,j,k,l]);
                       diff = abs(elecelecrepcalc[i,j,k,l]-elecelecreptest[i,j,k,l]);
                       @test diff <= tolerance;
	           end
               end #k=1:4 
           end
       end #i=1:4
        

    end#testset

end #run_elecelecrepul_test()


"""
    Test data

"""
function gettestdata()
   testdata = zeros(4,4,4,4);
    
   testdata[:,:,1,1]= [1.07657   0.57847   0.77610   0.53268;
   0.57847   0.58740   0.48970   0.55254;
   0.77610   0.48970   0.89630   0.51160;
   0.53268   0.55254   0.51160   0.55730];

   testdata[:,:,2,1] =[0.57847   0.32942   0.42735   0.30366;
   0.32942   0.36130   0.28576   0.34043;
   0.42735   0.28576   0.51160   0.30076;
   0.30366   0.34043   0.30076   0.34503];

   testdata[:,:,3,1] = [0.77610   0.42735   0.61361   0.40139;
   0.42735   0.44585   0.40139   0.43066;
   0.61361   0.40139   0.77610   0.42735;
   0.40139   0.43066   0.42735   0.44585];

   testdata[:,:,4,1] = [0.53268   0.30366   0.40139   0.28144;
   0.30366   0.33360   0.27019   0.31697;
   0.40139   0.27019   0.48970   0.28576;
   0.28144   0.31697   0.28576   0.32383];

   testdata[:,:,1,2] = [0.57847   0.32942   0.42735   0.30366;
   0.32942   0.36130   0.28576   0.34043;
   0.42735   0.28576   0.51160   0.30076;
   0.30366   0.34043   0.30076   0.34503];

   testdata[:,:,2,2] = [0.58740   0.36130   0.44585   0.33360;
   0.36130   0.45315   0.32383   0.42858;
   0.44585   0.32383   0.55730   0.34503;
   0.33360   0.42858   0.34503   0.43911];

   testdata[:,:,3,2] = [0.48970   0.28576   0.40139   0.27019;
   0.28576   0.32383   0.28144   0.31697;
   0.40139   0.28144   0.53268   0.30366;
   0.27019   0.31697   0.30366   0.33360];

   testdata[:,:,4,2] = [0.55254   0.34043   0.43066   0.31697;
   0.34043   0.42858   0.31697   0.41183;
   0.43066   0.31697   0.55254   0.34043;
   0.31697   0.41183   0.34043   0.42858];

   testdata[:,:,1,3] = [0.77610   0.42735   0.61361   0.40139;
   0.42735   0.44585   0.40139   0.43066;
   0.61361   0.40139   0.77610   0.42735;
   0.40139   0.43066   0.42735   0.44585];

   testdata[:,:,2,3] = [0.48970   0.28576   0.40139   0.27019;
   0.28576   0.32383   0.28144   0.31697;
   0.40139   0.28144   0.53268   0.30366;
   0.27019   0.31697   0.30366   0.33360];

   testdata[:,:,3,3] = [0.89630   0.51160   0.77610   0.48970;
   0.51160   0.55730   0.53268   0.55254;
   0.77610   0.53268   1.07657   0.57847;
   0.48970   0.55254   0.57847   0.58740];

   testdata[:,:,4,3] = [0.51160   0.30076   0.42735   0.28576;
   0.30076   0.34503   0.30366   0.34043;
   0.42735   0.30366   0.57847   0.32942;
   0.28576   0.34043   0.32942   0.36130];

   testdata[:,:,1,4] = [0.53268   0.30366   0.40139   0.28144;
   0.30366   0.33360   0.27019   0.31697;
   0.40139   0.27019   0.48970   0.28576;
   0.28144   0.31697   0.28576   0.32383];

   testdata[:,:,2,4] = [0.55254   0.34043   0.43066   0.31697;
   0.34043   0.42858   0.31697   0.41183;
   0.43066   0.31697   0.55254   0.34043;
   0.31697   0.41183   0.34043   0.42858];

   testdata[:,:,3,4] = [0.51160   0.30076   0.42735   0.28576;
   0.30076   0.34503   0.30366   0.34043;
   0.42735   0.30366   0.57847   0.32942;
   0.28576   0.34043   0.32942   0.36130];

   testdata[:,:,4,4] = [0.55730   0.34503   0.44585   0.32383;
   0.34503   0.43911   0.33360   0.42858;
   0.44585   0.33360   0.58740   0.36130;
   0.32383   0.42858   0.36130   0.45315];
 	   
    return testdata
end #gettestdata
