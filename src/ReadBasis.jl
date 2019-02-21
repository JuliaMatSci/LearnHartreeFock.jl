module ReadBasis

using BuildBasis,TypesBasis

export readgaussparams, retreadbasis




"""
	use: readgaussparams(basisfile,format="nwchem")
input:
    basisfile - name of file containing basis parameters.
    format (optional) - the format of the basis file.

returns: basisread - lines read from basis file

Description: Read the exponent and coefficients for a gaussian basis set. Returns an array of
strings that need to be later processed using TypesBasis.jl and Basis.jl.

Supports: mod NWChem EMSL format, i.e., but ignores species and orbital info by
commenting them out in the file with #.

""" function readgaussparams(filename;format="modnwchem")
            ele = geteledict();
            basisread = String[];
	    if format == "modnwchem"
	    	  file = open(filename,"r");
		  for line in readlines(file)
           	       if occursin("#",line)
               		   continue;
           	       elseif isempty(line)
               		    continue;
           	       else
			    push!(basisread,line);
           	       end
       		  end
		close(file);

                nbasisread = length(basisread);
                i = 1;
                nbasisfunc = 0;
                basisfunc = [];
                flag = false;
                while flag != true
                    line = basisread[i]
                    sline = split(line);
                    if haskey(ele,sline[1])
                        symbol = sline[1];
                        orbital = sline[2];
                        nprims = parse(Int,sline[3]);
                        nbasisfunc += 1;
                    else
                        error("Basis file needs start line having element oribtal number-primitives, example, H S 3")
                    end
                    params = [];
                    for ii=1:nprims
                        i += 1
                        push!(params,basisread[i]);
                    end
                    push!(basisfunc,params);                    
                    i += 1;
                    if i > nbasisread
                        flag = true;
                    end
                end # while flag
                
	        return nbasisfunc,basisfunc
             else
		error("Format $format is not supported yet.")
	     end 
end #readgaussparams 

"""
     use: retgaussbasis(myfile, format="nwchem", btype="STO-3G", atomicZ=1)
Input:
Returns: TypesBasis.BasisFunctions

Description: constructs TypeBasis.Gaussian.

""" function retreadbasis(filename;format="modnwchem",btype="6-31G")
    
    if isequal(btype,"6-31")
        #Parse file
        nbasisfunc,basisfuncstring = readgaussparams(filename,format=format);

        #Construct datatype
        basisfuncarry = Array{Gaussian,1}(undef,nbasisfunc);
        for (ibfs,bfs) in enumerate(basisfuncstring)
            nprims = length(bfs);
            #Unpack Array{String,1} to two Array{Float64,1}
            exponents = Array{Float64,1}(undef,nprims);
            coefficients = Array{Float64,1}(undef,nprims);
            for i=1:nprims
                strexpo,strcoef = split(bfs[i]);
                exponents[i] = parse(Float64,strexpo);
                coefficients[i] = parse(Float64,strcoef);
            end
            basisfuncarry[ibfs] = Gaussian(btype,nprims,coefficients,exponents)
        end #enumerate(basisfuncstring)        
        return BasisFunctions(nbasisfunc,basisfuncarry)
    else
        error("The basis set type $btype is not implemented")
    end 

end #retreadbasis

geteledict() = Dict("H" => 1, "He" => 2, "Li" => 3);

end #ReadBasis
