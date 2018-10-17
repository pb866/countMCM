#!/usr/local/bin/julia

# Load Julia packages
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
function readDB(filename::String="MCMv331species.db"; folder = "../DB")
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
    MCMspecies = strip.(replace.(lines[istart:iend], r"=[ ]*IGNORE[ ]*;" => ""))
  end

  return MCMspecies
end


"""
    translateSPC(species::AbstractString, old_name::String, new_name::String, spcDB::DataFrame=MCMdb, version::String="v3.3.1")

For a `species` name in the nomenclature `old_name` return a `String` with the `new_name`
by using the `spcDB` `DataFrame` with all translations and specifying the MCM `version`
as `v3.x`.
"""
function translateSPC(species::AbstractString, old_name::String, new_name::String,
  spcDB::DataFrame=MCMdb, version::String="v3.3.1")
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


# MCMv3.3.1
MCMdb = readDB()
MCMspecies = readKPPspc()
MCMgecko = translateSPC.(MCMspecies, "MCMname", "GECKO-A")
