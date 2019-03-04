module CoulombExchange

#LearnHatreeFock.jl Modules
using TypesBasis,TypesParticles


export getcoulexchange

@doc raw"""
function getcoulexchange()

"""
function getcoulexchange(density::Array,elecelecrepul::Array)

    nbasisfunc = size(density)[1];
    coulexchange = zeros(Float64,nbasisfunc,nbasisfunc);

    for l=1:nbasisfunc
        for m=1:nbasisfunc
            for n=1:nbasisfunc
                for o=1:nbasisfunc
                    diff = 2.00e0*elecelecrepul[l,m,n,o]-elecelecrepul[l,n,m,o];
                    coulexchange[l,m] += density[n,o]*(diff);
                end
            end #n=1:nbasisfunc
        end
    end #l=1:nbasisfunc
    
    return coulexchange
end #getcoulexchange


end #CoulombExchange
