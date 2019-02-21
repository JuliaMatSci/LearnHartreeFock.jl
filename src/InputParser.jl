module InputParser

using TypesParticles, TypesBasis
using ReadBasis


export getcalcsetup,parseinputfile

"""
    Use: getcalcsetup(filename::String)

Input: 
    filename - input file name and location

Return:
    system - TypesParticles.SystemOfAtoms for calculation
    basis - TypesBasis for calculation

Description: this function is called by the main driver program HatreeFock.jl to
setup the calculation.

TODO: Unit test functions needed

""" function getcalcsetup(filename::String)
    results = parseinputfile(filename);

    system = SystemOfAtoms(results["natoms"],
                           results["positions"],
                           results["cell"]);
    btype = results["basistype"];
    if btype == "6-31"
        basis = retreadbasis(results["basisfile"],btype=btype);
    else
        error("Basis type $btype is not implemented.")
    end
    
    return (system,basis)
end #getcalcsetup

"""
        Use: parseinputfile(filename)
Inputs: 
        filename - input file name and location
Returns:
        parserresults - dictionary of keywords => function evaluations

Description: This function parses the input file and then calls the respective keyword
get function to parse the line. The function returns a dictionary whos keys are
the input keywords and values are the results from the function call.


Input syntax: keyword list-of-commands
                   & - line read continuation character
                   # - line read comment character, included in function call return.

Notes: All errors associated with keyword lines are produced in get_xxxx calls

TODOs: unit test to ensure all gotchas

""" function parseinputfile(filename::String)

    if ! isfile(filename)
        error("The file $filename does not exist")
    end

    parseresults = Dict();
    infile = open(filename,"r")
    keywords = getreqkeywords();
    comments=Array{String,1}();


    while !eof(infile)
        line = readline(infile);
        if (occursin("#",line) || isempty(line))
            push!(comments,line);
        else
            for k in keywords
            #Construct function name to call based on keyword
            funcname = Symbol(:get,k);
    
            if occursin(k,line)
                linestring = split(line,"$k ")[end];

                # Read keyword multiple lines using \& symbol
                while occursin(r" * &\s?$",linestring)
                    nextline=readline(infile);
                    linestring = linestring*nextline;
                end

                #CHECK - do I want the line below, requires complex string parsing if so.
                linestring=replace(linestring, "&" => "");
    
                #Evaluate and insert results into dictionary
                parseresults[k] = eval(funcname)(linestring);
            end
            end #keywords
        end
    end #eof

   #Need to loop over any required keywords that were not found in file.
   #TODO - I assume a more clever approach could be used  here.
    for k in keywords
        try
            parseresults[k]
        catch
            funcname = Symbol(:get,k);
            parseresults[k] = eval(funcname)("");
        end
    end
    
    parseresults["comments"] = comments;
    close(infile)
    
    return parseresults

end #parseinputfile

"""
        Use: get_cell(line)

Inputs: 
    line - String that can be unpacked to a (3,1) Array
    default - Specify if keyword is found was found in read file.
Returns: 
    default = [10.0; 10.0; 10.0] or
    3-element Array{Float64,1}

Description: parse simulation x, y, and z lengths in Angstroms. Only
supports cubic cells.


""" function getcell(line::String)
    splitline = split(line);

    
    if isempty(line)
        println("NOTICE: default cell 10.0 10.0 10.0 is used!")
        return [10.00e0 0.00e0 0.00e0;
                0.00e0 10.00e0 0.00e0;
                0.00e0 0.00e0 10.00e0]
     elseif length(splitline) != 3
            error("The length of the parsed cell string does not correspond to 3-element array.")
    else
        #Only supports cubic cell
        diag = parse.(Float64,split(line))
        cell = zeros(3,3);
        cell[1,1],cell[2,2],cell[3,3] = diag[1], diag[2], diag[3];
        return cell
    end 
end #get_cell


"""
        Use: get_positions(lines::String
Inputs: lines - is the atomic number and positions for 1 or more atoms
Returns: atoms - Array of data type Atom 


""" function getpositions(lines::String)

    nslots = 4;
    atoms = Array{Atom,1}();
    if isempty(lines)
        println("NOTICE: using default position 5.0 5.0 5.0 for 1 hydrogen atom!")
        return push!(atoms,Atom(5.00e0,5.00e0,5.00e0,1))
    else
        splitlines = split(lines);
        natoms = Int(length(splitlines)/nslots);
        nentries = reshape(splitlines,(nslots,natoms));
        for n=1:natoms
            Z = parse(Int,nentries[1,n]);
            x,y,z = parse.(Float64,nentries[2:4,n]);
            push!(atoms,Atom(x,y,z,Z));
        end #natoms
        return atoms
    end 
    
end #get_positions


"""
        Use: get_natoms(line::String)
Inputs: line - number of atoms 
Returns: Integer data type of line

""" function getnatoms(line::String)
    if isempty(line)
        println("NOTICE: using 1 atom as default!")
        return 1
    elseif length(replace(line," " => "")) != 1
        println(line)
        error("The length of the parsed natoms string does not correspond to a single value.")
    else    
        return parse(Int,line)
    end 
    
end #get_natoms




"""
returns allowed keywords in input file.
""" function getreqkeywords()
    keywords = Set{String}(["basisfile",
                            "basistype",
                            "cell",
                            "natoms",
                            "positions"
                            ])
end #getkeywords


#Inline functions 
getbasisfile(line::String)=line;
getbasistype(line::String)=replace(line," " => "");



end #InputParser
