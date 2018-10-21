#!/usr/local/bin/julia

# Load Julia packages
using Pkg
Pkg.activate(normpath(joinpath(Base.source_dir(),"..")))
using DataFrames
#=
Solutions to precomile errors:
• Build package "CodecZlib" (run `build "CodecZlib"` in the package manager)
=#


"""
    readDB(filename::String="MCMspecies.db"; folder::String="DB")

Read database file with `filename` in default or specified `folder` and return
a DataFrame with columns for MCM names, SMILES, InChIs, GECKO-A names, and molar masses.

Column headers are defined in the first line of the database file.
Columns are separated by an ampersand sign (`&`), use `&` also at the end of each line
as last character. Missing species names can be put in curly braclets (`{...}`) and will
be ignored in the database.
"""
function readDB(filename::String="MCMspecies.db"; folder::String="DB")
  MCMdb = DataFrame()
  open(joinpath(pwd(), folder, filename)) do f
    lines = readlines(f)
    species_data = split.(lines, "&")
    for i = 1:length(species_data[1])-1
      MCMdb[Symbol(species_data[1][i])] = [dat[i] for dat in species_data[2:end]]
    end
  end

  return MCMdb
end


"""
    readKPPspc(KPPfile::String; folder::String = "KPPfiles")

From an MCM `KPPfile` in a default or specified `folder` read all species names
(in MCM nomenclature) and return as array of strings.
"""
function readKPPspc(KPPfile::String; folder::String = "KPPfiles")
  MCMspecies = String[]
  open(joinpath(pwd(), folder, KPPfile)) do f
    lines = readlines(f)
    idx   = findall(occursin.("IGNORE", lines))
    MCMspecies = strip.([replace(l, r"=[ ]*IGNORE[ ]*;" => "") for l in lines[idx]])
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
    readRO2(KPPfile::String="MCMv3.3.1.kpp"; folder::String="KPPfiles")

Retrieve all RO2 species (MCM names) from a default or specified `folder`/`KPPfile`
and return the names as a Vector of strings.
"""
function readRO2(KPPfile::String; folder::String="KPPfiles")
  # Initialise RO2 list
  RO2 = String[]
  open(joinpath(pwd(), folder, KPPfile), "r") do f
    # Read KPP file and find section with RO2 summation
    lines = readlines(f)
    RO2index = findall(occursin.("C(", lines))
    RO2lines = split.([replace(l, "&" => "") for l in lines[RO2index]], "+")
    # Add RO2 species names to a list with out the KPP syntax for concentrations
    for line in RO2lines, spc in line
      spc = replace(spc, "C(ind_" => "")
      spc = strip.(replace(spc, ")" => ""))
      if spc ≠ ""  push!(RO2, spc)  end
    end
  end

  return RO2
end


"""
    findRO2(MCMspecies::Union{Vector{String}, Vector{SubString{String}}}, old_name::String="MCMname", new_name::String="GECKO-A"; spcDB::DataFrame=MCMdb, version::String="v3.3.1")

From a list of `MCMspecies`, return all RO2 as vector of strings temporarily translating
from `old_name` (`"MCMname"`) to `new_name` (`"GECKO-A"`) using the translation database
`spcDB` for the specified MCM `version` 3.x.
"""
function findRO2(MCMspecies::Union{Vector{String}, Vector{SubString{String}}},
  old_name::String="MCMname", new_name::String="GECKO-A";
  spcDB::DataFrame=MCMdb, version::String="v3.3.1")
  RO2 = String[]
  for ro2 in MCMspecies
    ro2gecko = translateSPC(ro2, old_name, new_name, spcDB = spcDB, version = version)
    if occursin("(OO.)", ro2gecko) && !occursin(".(OO.)", ro2gecko)
      push!(RO2, ro2)
    end
  end

  return RO2
end


"""
    RO2conflicts(RO2searchDB::Vector{String}, RO2SPClist::Vector{String})

Find all conflicts for RO2 in the KPP sum or the mechanism from a list of RO2 (`RO2SPClist`)
compared to a list (`RO2RO2searchDB`).
"""
function RO2conflicts(RO2searchDB::Vector{String}, RO2SPClist::Vector{String})
  RO2returnDB = String[]
  for ro2 in RO2SPClist
    if all(RO2searchDB .≠ ro2)
      println("RO₂ conflict: $ro2")
      push!(RO2returnDB, ro2)
    end
  end

  return RO2returnDB
end



### MCMv3.3.1 ###
# Read databases and KPP files to find and analyse RO2 in the MCM mechanism and the RO2 sum
MCMdb = readDB()
MCMv33species = readKPPspc("MCMv3.3.1.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm33 = findRO2(MCMv33species)
RO2sum33 = readRO2("MCMv3.3.1.kpp")

# Find missing RO2 and RO2 conflicts
RO2v33missing = RO2conflicts(RO2sum33, RO2mcm33)
RO2v33conflicts = RO2conflicts(RO2mcm33, RO2sum33)


### MCMv3.2 ###
# Find current species in mechanism
MCMv32species = readKPPspc("MCMv3.2.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm32 = findRO2(MCMv32species, version = "MCMv3.2")
RO2sum32 = readRO2("MCMv3.2.kpp")

# Find missing RO2 and RO2 conflicts
RO2v32missing = RO2conflicts(RO2sum32, RO2mcm32)
RO2v32conflicts = RO2conflicts(RO2mcm32, RO2sum32)


### MCMv3.1 ###
# Find current species in mechanism
MCMv31species = readKPPspc("MCMv3.1.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm31 = findRO2(MCMv31species, version = "MCMv3.1")
RO2sum31 = readRO2("MCMv3.1.kpp")

# Find missing RO2 and RO2 conflicts
RO2v31missing = RO2conflicts(RO2sum31, RO2mcm31)
RO2v31conflicts = RO2conflicts(RO2mcm31, RO2sum31)


# Print conflicts to file
open("RO2conflicts.dat", "w") do f
  if !isempty(RO2v33missing)
    println(f, "The following RO2 are missing in the summation of MCMv3.3.1:")
    println(f, join(RO2v33missing, ", "),"\n")
  end
  if !isempty(RO2v33conflicts)
    println(f, "The following RO2 should not be in the summation of MCMv3.3.1:")
    println(f, join(RO2v33conflicts, ", "),"\n")
  end
  if !isempty(RO2v32missing)
    println(f, "The following RO2 are missing in the summation of MCMv3.2:")
    println(f, join(RO2v32missing, ", "),"\n")
  end
  if !isempty(RO2v32conflicts)
    println(f, "The following RO2 should not be in the summation of MCMv3.2:")
    println(f, join(RO2v32conflicts, ", "),"\n")
  end
  if !isempty(RO2v31missing)
    println(f, "The following RO2 are missing in the summation of MCMv3.1:")
    println(f, join(RO2v31missing, ", "),"\n")
  end
  if !isempty(RO2v31conflicts)
    println(f, "The following RO2 should not be in the summation of MCMv3.1:")
    println(f, join(RO2v31conflicts, ", "),"\n")
  end
end

println("done.")
