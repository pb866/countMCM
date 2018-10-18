#!/usr/local/bin/julia

# Load Julia packages
using Pkg
Pkg.activate(".")
using DataFrames
#=
Solutions to precomile errors:
• Build package "CodecZlib" (run `build "CodecZlib"` in the package manager)
=#


"""
    readDB(filename::String="MCMv331species.db"; folder = "../DB")

Read database file with `filename` in default or specified `folder` and return
a DataFrame with columns for MCM names, SMILES, InChIs, GECKO-A names, and molar masses.

Column headers are defined in the first line of the database file.
Columns are separated by an ampersand sign (`&`), use `&` also at the end of each line
as last character. Missing species names can be put in curly braclets (`{...}`) and will
be ignored in the database.
"""
function readDB(filename::String="MCMspecies.db"; folder = "../DB")
  MCMdb = DataFrame()
  open(joinpath(Base.source_dir(), folder, filename)) do f
    lines = readlines(f)
    species_data = split.(lines, "&")
    for i = 1:length(species_data[1])-1
      MCMdb[Symbol(species_data[1][i])] = [dat[i] for dat in species_data[2:end]]
    end
  end

  return MCMdb
end


"""
    readKPPspc(KPPfile::String="MCMv3.3.1.kpp"; folder = "../KPPfiles")

From an MCM `KPPfile` in a default or specified `folder` read all species names
(in MCM nomenclature) and return as array of strings.
"""
function readKPPspc(KPPfile::String="MCMv3.3.1.kpp"; folder = "../KPPfiles")
  MCMspecies = String[]
  open(joinpath(Base.source_dir(), folder, KPPfile)) do f
    lines = readlines(f)
    istart = findfirst(occursin.("#DEFVAR", lines)) + 1
    iend   = findlast(occursin.("IGNORE", lines))
    MCMspecies = strip.([replace(l, r"=[ ]*IGNORE[ ]*;" => "") for l in lines[istart:iend]])
  end

  return MCMspecies
end


"""
    translateSPC(species::AbstractString, old_name::String, new_name::String, spcDB::DataFrame=MCMdb, version::String="v3.3.1")

For a `species` name in the nomenclature `old_name` return a `String` with the `new_name`
by using the `spcDB` `DataFrame` with all translations and specifying the MCM `version`
as `v3.x`.
"""
function translateSPC(species::AbstractString, old_name::String, new_name::String;
  spcDB::DataFrame=MCMdb, version::String="v3.3.1")
  if species == ""  return "DUMMY"  end # Return empty string as DUMMY
  if version == "v3.3.1"
    i = findfirst(spcDB[Symbol(old_name)] .== species)
  else
    i = findlast(spcDB[Symbol(old_name)] .== species)
  end
  return_species = nothing
  if i == nothing
    println("$species not found! Translation skipped.")
    println("Add species to translation database.")
  else
    return_species = spcDB[Symbol(new_name)][i]
  end

  return return_species
end


"""
    readRO2(KPPfile::String="MCMv3.3.1.kpp"; folder::String="../KPPfiles")

Retrieve all RO2 species (MCM names) from a default or specified `folder`/`KPPfile`
and return the names as a Vector of strings.
"""
function readRO2(KPPfile::String="MCMv3.3.1.kpp"; folder::String="../KPPfiles")
  # Initialise RO2 list
  RO2 = String[]
  open(joinpath(Base.source_dir(), folder, KPPfile), "r") do f
    # Read KPP file and find section with RO2 summation
    lines = readlines(f)
    RO2index = findall(occursin.("C(", lines))
    RO2lines = split.([replace(l, "&" => "") for l in lines[RO2index]], "+")
    # Add RO2 species names to a list with out the KPP syntax for concentrations
    for line in RO2lines, spc in line
      println(spc)
      spc = replace(spc, "C(ind_" => "")
      spc = strip.(replace(spc, ")" => ""))
      if spc ≠ ""  push!(RO2, spc)  end
    end
  end

  return RO2
end


### MCMv3.3.1 ###
# Read databases and KPP files to find and analyse RO2 in the MCM mechanism and the RO2 sum
MCMdb = readDB()
MCMv33species = readKPPspc()
MCMv33gecko = translateSPC.(MCMv33species, "MCMname", "GECKO-A")
RO2 = findall(occursin.("(OO.)", MCMv33gecko))

# Define RO2 in the RO2 sum and the mechanism
RO2mcm33 = MCMv33species[RO2]
RO2sum33 = readRO2()


### MCMv3.2 ###
# Read databases and KPP files to find and analyse RO2 in the MCM mechanism and the RO2 sum
MCMv32species = readKPPspc("MCMv3.2.kpp")
MCMv32gecko = translateSPC.(MCMv32species, "MCMname", "GECKO-A", version = "v3.2")
RO2 = findall(occursin.("(OO.)", MCMv32gecko))

# Define RO2 in the RO2 sum and the mechanism
RO2mcm32 = MCMv32species[RO2]
RO2sum32 = readRO2("MCMv3.2.kpp")

### MCMv3.1 ###
# Read databases and KPP files to find and analyse RO2 in the MCM mechanism and the RO2 sum
MCMv31species = readKPPspc("MCMv3.1.kpp")
MCMv31gecko = translateSPC.(MCMv31species, "MCMname", "GECKO-A", version = "v3.1")
RO2 = findall(occursin.("(OO.)", MCMv31gecko))

# Define RO2 in the RO2 sum and the mechanism
RO2mcm31 = MCMv31species[RO2]
RO2sum31 = readRO2("MCMv3.1.kpp")

println("done.")
