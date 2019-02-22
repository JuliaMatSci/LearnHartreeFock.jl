module UtilityFunctions

export flattenbasisfunc

using TypesParticles, TypesBasis


"""
TODO write description

""" function flattenbasisfunc(numelec::Int,numbasisfunc::Int,
                              basisfunc::Array{GaussOrbitals})
    overlapsize = numelec*numbasisfunc;
    flatbasisfunc = Array{GaussOrbitals,1}(undef,overlapsize);
    ei = 1;
    for e=1:numelec
        tmpbf = basisfunc[e,:]
        for i=1:numbasisfunc
            flatbasisfunc[ei] = tmpbf[i];
            ei += 1;
        end
    end
    return flatbasisfunc
end #flattenbasisfunc


end #module
