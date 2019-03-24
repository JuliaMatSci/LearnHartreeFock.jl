export LEARNHATREEFOCK_PATH="/home/stefanb/GoogleDrive/COMPUTATION/LearnHatreeFock.jl/src"

echo "Distance[Bohr] Energy[Ha]" > He2.energy_distance.dat

distances="0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0 2.2 2.4 2.6 2.8 3.0 3.4 3.8 4.4"
for pos in ${distances}
do
    cat  <<EOF > He2.in

    #Hydrogen molecule

    #Angstroms
    cell 10.0 10.0 10.0

    natoms 2

    #Atoms Z ,x , y, z in Angstroms
    positions 2 0.00 0.00 0.00 &
    2 ${pos} 0.00 0.00 

    #Basis definition
    basistype 6-31
    basisfile ../../examples/basissets/He.6-31G.mod
    
EOF

    julia ../../src/HatreeFock.jl He2.in > He2.out
    energy=`grep 'Total Minimum Energy: *' He2.out | awk '{print $4}'`
    echo ${pos} ${energy} >> He2.energy_distance.dat
done

