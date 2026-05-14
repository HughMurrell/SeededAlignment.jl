using Pkg
Pkg.activate("../.")
using FASTX, BioSequences, StatsBase, DataFrames, CSV
using SeededAlignment


# test script starts here
# to try it out type:
# julia FunctionalFilter.jl consensusC_env_nt.fasta CAP304_420.fasta 0.8

if length(ARGS)!=3
    println("usage: julia FunctionalFilter.jl in_file, ref_file min_match")
    exit()
end
ref_file=ARGS[1]
in_file=ARGS[2]
min_match=parse(Float64,ARGS[3])
@show min_match
hk = filter_and_align(ref_file,in_file,"functionals.fasta","nonfunctionals.fasta", match_thresh=min_match)
@show hk

