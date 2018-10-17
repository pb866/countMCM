using DataFrames
#=
Solutions to precomile errors:
• Build package "CodecZlib" (run `build "CodecZlib"` in the package manager)
=#

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

# MCMv3.3.1
MCMdb = read()
MCMspecies = readKPPspc()
