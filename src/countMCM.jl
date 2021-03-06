#!/usr/local/bin/julia

# Load Julia packages
using Pkg
Pkg.activate(normpath(joinpath(Base.source_dir(),"..")))
using DataFrames
using filehandling
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
    idx   = findall(MCMspecies.=="")
    deleteat!(MCMspecies, idx)
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
    readRO2(KPPfile::String; folder::String="KPPfiles")

Retrieve all RO2 species (MCM names) from a `KPPfile` in a default or specified `folder`
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
    conflicts(RO2searchDB::Union{Vector{String}, Vector{SubString{String}}}, RO2SPClist::Union{Vector{String}, Vector{SubString{String}}})

Find all conflicts for RO2 in the KPP sum or the mechanism from a list of RO2 (`RO2SPClist`)
compared to a list (`RO2RO2searchDB`).
"""
function conflicts(searchDB::Union{Vector{String}, Vector{SubString{String}}},
  SPClist::Union{Vector{String}, Vector{SubString{String}}})
  returnDB = String[]
  for spc in SPClist
    if all(searchDB .≠ spc)
      println("Species conflict: $spc")
      push!(returnDB, spc)
    end
  end

  return returnDB
end


"""
checkFAC(FACfile::String, MCMspecies::Union{Vector{String}, Vector{SubString{String}}}, RO2mcm::Union{Vector{String}, Vector{SubString{String}}}; folder::String="FACfiles")

From a `FACfile` in a default or specified `folder`, read all species and RO2 definitions
and compare them to the arrays of `MCMspecies` and `RO2mcm` with the KPP file definitions.

Returns `String` arrays with FAC species and RO2 names and a tuple with booleans,
whether differences where found for species and RO2, respectively.
"""
function checkFAC(FACfile::String, MCMspecies::Union{Vector{String}, Vector{SubString{String}}},
  RO2mcm::Union{Vector{String}, Vector{SubString{String}}}; folder::String="FACfiles")
  # Initialise
  FACspecies = String[]; RO2fac = String[]
  SPCflag = true; RO2flag = true
  open(joinpath(pwd(),folder,FACfile)) do f
    # Read FAC file
    lines = readlines(f)
    # Find species definitions
    varstart = findfirst(lines.=="VARIABLE") + 1
    varend = findnext(occursin.(";", lines), varstart)
    lines[varend] = replace(lines[varend], ";" => "")
    FACspc = split.(lines[varstart:varend])
    for line in FACspc, spc in line
      push!(FACspecies, strip(spc))
    end
    # Compare to KPP file and warn, if different
    SPCflag = sort(MCMspecies) ≠ sort(FACspecies)
    if SPCflag
      println("WARNING! Different species lists for KPP and FAC file!")
    end
    # Find RO2 definitions
    ro2start = findfirst([occursin(r"RO2[ ]*=", l) for l in lines])
    ro2end = findnext(occursin.(";", lines), ro2start)
    lines[ro2start] = replace(lines[ro2start], r"RO2[ ]*=" => "")
    lines[ro2end] = replace(lines[ro2end], ";" => "")
    for line in lines[ro2start:ro2end]
      ro2line = split(line,"+")
      for ro2 in ro2line
        if strip(ro2) ≠ ""
          push!(RO2fac, strip(ro2))
        end
      end
    end
    # Compare to KPP file and warn, if different
    RO2flag = sort(RO2fac) ≠ sort(RO2mcm)
    if RO2flag
      println("WARNING! Different RO2 lists for KPP and FAC file!")
    end
  end

  return FACspecies, RO2fac, (SPCflag, RO2flag)
end


"""
    FACvsKPP(flag::Tuple{Bool, Bool}, MCMspecies::Union{Vector{String}, Vector{SubString{String}}},
      FACspecies::Union{Vector{String}, Vector{SubString{String}}},
      MCMro2::Union{Vector{String}, Vector{SubString{String}}},
      FACro2::Union{Vector{String}, Vector{SubString{String}}},
      version::String)

If the respective `flag` is true, compare the respective MCM and FAC species and/or RO2
and print missing species in either file to a warning output file labelled with the
specified `version` of the MCM mechanism.
"""
function FACvsKPP(flag::Tuple{Bool, Bool}, MCMspecies::Union{Vector{String}, Vector{SubString{String}}},
  FACspecies::Union{Vector{String}, Vector{SubString{String}}},
  MCMro2::Union{Vector{String}, Vector{SubString{String}}},
  FACro2::Union{Vector{String}, Vector{SubString{String}}},
  version::String)

  if flag[1]
    SPCmissing  = conflicts(FACspecies, MCMspecies)
    SPCconflict = conflicts(MCMspecies, FACspecies)
    open("conflicts_species_MCM$version.dat", "w") do f
      if !isempty(SPCmissing)
        println(f, "The following KPP species definitions are missing in the FAC file:")
        println(f, join(SPCmissing, ", "), '\n')
      end
      if !isempty(SPCconflict)
        println(f, "The following FAC species definitions are missing in the KPP file:")
        println(f, join(SPCconflict, ", "), '\n')
      end
    end
  end
  if flag[2]
    RO2missing  = conflicts(FACro2, MCMro2)
    RO2conflict = conflicts(MCMro2, FACro2)
    open("conflicts_RO2_MCM$version.dat", "w") do f
      if !isempty(RO2missing)
        println(f, "The following KPP species definitions are missing in the FAC file:")
        println(f, join(RO2missing, ", "), '\n')
      end
      if !isempty(RO2conflict)
        println(f, "The following FAC species definitions are missing in the KPP file:")
        println(f, join(RO2conflict, ", "), '\n')
      end
    end
  end
end


### MCMv3.3.1 ###
println("Checking MCMv3.3.1...")
# Read databases and KPP files to find and analyse RO2 in the MCM mechanism and the RO2 sum
MCMdb = readDB()
MCMv33species = readKPPspc("MCMv3.3.1.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm33 = findRO2(MCMv33species)
RO2sum33 = readRO2("MCMv3.3.1.kpp")

# Find missing RO2 and RO2 conflicts
RO2v33missing = conflicts(RO2sum33, RO2mcm33)
RO2v33conflicts = conflicts(RO2mcm33, RO2sum33)
# Compare KPP and FAC files
FACv33species, FACv33ro2, v33flag = checkFAC("MCMv3.3.1.fac", MCMv33species, RO2sum33)
FACvsKPP(v33flag, MCMv33species, FACv33species, FACv33ro2, RO2sum33, "v3.3.1")


### MCMv3.2 ###
println("Checking MCMv3.2...")
# Find current species in mechanism
MCMv32species = readKPPspc("MCMv3.2.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm32 = findRO2(MCMv32species, version = "MCMv3.2")
RO2sum32 = readRO2("MCMv3.2.kpp")

# Find missing RO2 and RO2 conflicts
RO2v32missing = conflicts(RO2sum32, RO2mcm32)
RO2v32conflicts = conflicts(RO2mcm32, RO2sum32)
# Compare KPP and FAC files
FACv32species, FACv32ro2, v32flag = checkFAC("MCMv3.2.fac", MCMv32species, RO2sum32)
intersect(MCMv32species, FACv32species)
FACvsKPP(v32flag, MCMv32species, FACv32species, FACv32ro2, RO2sum32, "v3.2")


### MCMv3.1 ###
println("Checking MCMv3.1...")
# Find current species in mechanism
MCMv31species = readKPPspc("MCMv3.1.kpp")

# Define RO2 in the RO2 sum and the mechanism
RO2mcm31 = findRO2(MCMv31species, version = "MCMv3.1")
RO2sum31 = readRO2("MCMv3.1.kpp")
RO2sum31Leeds = readRO2("MCMv3.1_Leeds.kpp")

# Find missing RO2 and RO2 conflicts
RO2v31missing = conflicts(RO2sum31, RO2mcm31)
RO2v31conflicts = conflicts(RO2mcm31, RO2sum31)
RO2v31missingLeeds = conflicts(RO2sum31Leeds, RO2mcm31)
RO2v31conflictsLeeds = conflicts(RO2mcm31, RO2sum31Leeds)
# Compare KPP and FAC files
FACv31species, FACv31ro2, v31flag =
  checkFAC("MCMv3.1.fac", MCMv31species, RO2sum31)
FACvsKPP(v31flag, MCMv31species, FACv31species, FACv31ro2, RO2sum31, "v3.1 (Legacy)")
FACv31speciesLeeds, FACv31ro2Leeds, v31flagLeeds =
  checkFAC("MCMv3.1.fac", MCMv31species, RO2sum31Leeds)
FACvsKPP(v31flagLeeds, MCMv31species, FACv31species, FACv31ro2, RO2sum31Leeds, "v3.1 (Leeds)")

println("Print results to output files \'RO2conflicts.dat\`.")
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
  if !isempty(RO2v31missingLeeds)
    println(f, "The following RO2 are missing in the summation of MCMv3.1:")
    println(f, join(RO2v31missingLeeds, ", "),"\n")
  end
  if !isempty(RO2v31conflictsLeeds)
    println(f, "The following RO2 should not be in the summation of MCMv3.1:")
    println(f, join(RO2v31conflictsLeeds, ", "),"\n")
  end
end


# Compare to AtChem
AtChem = readfile("FACfiles/peroxy-radicals_v3.1")

# Find missing RO2 and RO2 conflicts
AtChemLeeds = conflicts(RO2sum31Leeds, AtChem)
LeedsAtChem = conflicts(AtChem, RO2sum31Leeds)

println("There are $(length(AtChemLeeds)) additional RO₂ in AtChem compared to ",
  "the MCMv3.1 (Leeds).")
println("There are $(length(LeedsAtChem)) additional RO₂ in the MCMv3.1 (Leeds) ",
  "compared to AtChem.")


println("done.")
