function degap(s::String)
    return replace(s,"-"=>"")
end

function degap(s::LongDNA{4})
    return filter!(!isgap,s)
end

function remove_ambigs(s::String)
    return replace(s,"N"=>"")
end
     
function remove_ambigs(s::LongDNA{4})
    return filter!(!isambiguous,s)
end
     
function longest_open_reading_frame(cons)
    cons=degap(cons)
    cons_frames = [ cons[1:3*div(length(cons),3)],
                    cons[2:3*div(length(cons),3)-2],
                    cons[3:3*div(length(cons),3)-1] ]
    # then translate them
    cons_aa_frames = (x->String(BioSequences.translate(LongDNA{4}(degap(x))))).(cons_frames)
    cons_aa = cons_aa_frames[1]
        
    # find longest coding region in all frames
    best_reading_frame = 1
    start=1
    stop=length(cons_aa)
    max_match_length=0
    for rf in 1:3
        cons_aa = cons_aa_frames[rf]
        m=match(r"M[^\*]*\*",cons_aa)
        while ! isnothing(m)
            if length(m.match)>max_match_length
                best_reading_frame=rf
                start=m.offset
                stop=start-1+length(m.match)
                max_match_length=length(m.match)
            end
            m=match(r"M[^\*]*\*",cons_aa,m.offset+1)
        end
    end
    # cons_aa = cons_aa_frames[best_reading_frame]
    # cons_aa_trim = cons_aa[start:stop]
    cons_trim = cons[((start-1)*3+best_reading_frame):((stop)*3+best_reading_frame-1)]
    return cons_trim
end

function filter_and_align(ref_file, query_file, functionals_file, nonfunctionals_file;
            match_thresh=0.7)
    hk=DataFrame(sample=String[], sequences=Int[], functional=Int[], nonfunctional=Int[], percentLost=Int[],
        ambiguous=Int[],frameshift=Int[], lateStartCodon=Int[], earlyStopCodon=Int[],
        badMatch=Int[], matchThresh=Float64[])
    score_params = ScoringScheme(edge_ext_begin = true, edge_ext_end = true )
    ref_nams, ref_seqs = read_fasta(ref_file)
    ref_nam=ref_nams[1]
    ref_seq=ref_seqs[1]
    nams, seqs = read_fasta(query_file)
    seqs=(degap).(seqs)
    start_count = length(seqs)
    reject_seqs=[]
    reject_nams=[]
    keeps=(q->all(x -> x in (DNA_A, DNA_T, DNA_C, DNA_G), q)).(seqs)
    ambig_count=length(seqs)-sum(keeps)
    @show query_file, ambig_count
    reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*" ambiguousSymbols-reject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    seqs=seqs[keeps]
    for i in 1:length(seqs)
        seq_trim = longest_open_reading_frame(seqs[i])
        pw_align = seed_chain_align(ref=ref_seq, query=seq_trim, scoring=score_params, verbose=false)
        startTrim=findfirst('-'.!=(collect(string(pw_align[1]))))
        stopTrim=findlast('-'.!=(collect(string(pw_align[1]))))
        seqs[i]=filter!(!isgap, pw_align[2][startTrim:stopTrim])
    end
    keeps=(x->count(==(DNA_N),collect(x))==0).(seqs)
    orf_reject_count=length(seqs)-sum(keeps)
    reject_seqs=vcat(reject_seqs,seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*" frameshift-reject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    seqs=seqs[keeps]
    nams=vcat([ref_nam],nams)
    ali_seqs = msa_codon_align(ref_seq, seqs, scoring=score_params, verbose=false)
    startTrim=findfirst('-'.!=(collect(string(ali_seqs[1]))))
    stopTrim=findlast('-'.!=(collect(string(ali_seqs[1]))))
    trim_ali_seqs=(x->x[startTrim:stopTrim]).(ali_seqs)
    keeps=(x->(x[1:3]==dna"ATG")).(trim_ali_seqs)
    no_start_codon_count=length(trim_ali_seqs)-sum(keeps)
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*" lateStart-reject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    keeps=(x->(x[end-2:end]!=dna"---")).(trim_ali_seqs)
    no_stop_codon_count=length(trim_ali_seqs)-sum(keeps)
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*" earlyStop-reject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    match_ratios=(x->sum(collect(trim_ali_seqs[1]).==(collect(x)))/length(x)).(trim_ali_seqs)
    @show query_file, minimum(match_ratios), maximum(match_ratios), mean(match_ratios)
    keeps=match_ratios.>=match_thresh
    bad_match_count=sum((!).(keeps))
    reject_seqs=vcat(reject_seqs,trim_ali_seqs[(!).(keeps)])
    reject_nams=vcat(reject_nams,(x->x*" badMatch-reject").(nams[(!).(keeps)]))
    nams=nams[keeps]
    trim_ali_seqs=trim_ali_seqs[keeps]
    end_count = length(trim_ali_seqs)
    seq_loss = floor( Int, 100 * (start_count - end_count + 1) / start_count )
    @show query_file, seq_loss
    if length(trim_ali_seqs) > 1
        trim_ali_seqs = msa_codon_align(ref_seq, remove_ambigs.(degap.(trim_ali_seqs[2:end])),
                                            scoring=score_params, verbose=false)
    end
    write_fasta(functionals_file,trim_ali_seqs,seq_names=nams)
    if length(reject_seqs) > 0
        # write_fasta(in_file*"_functionalrejects.fasta",degap.(reject_seqs),seq_names=reject_nams)
        write_fasta(nonfunctionals_file,LongDNA{4}.(reject_seqs),seq_names=reject_nams)
    end
    hk_rec=[basename(query_file),start_count,end_count-1,start_count-end_count+1,seq_loss,
        ambig_count,orf_reject_count,no_start_codon_count,no_stop_codon_count,bad_match_count,match_thresh]
    push!(hk,hk_rec)
    return hk
end

"""
# test script starts here
# to try it out type:
# julia FunctionalFilter.jl CAP304_420.fasta consensusC_env_nt.fasta

if length(ARGS)!=3
    println("usage: julia FunctionalFilter.jl in_file, ref_file min_match")
    exit()
end
in_file=ARGS[1]
ref_file=ARGS[2]
min_match=parse(Float64,ARGS[3])
@show min_match
hk=DataFrame(sample=[], sequences=[], functional=[], non_functional=[],
                    frameshift_error=[], late_start_codon=[], early_stop_codon=[], bad_match=[])
pairwise_filter(in_file,ref_file,min_match_ratio=min_match,hk=hk)
@show hk
"""
