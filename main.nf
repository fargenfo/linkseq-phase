#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.bam_paths = null
params.vcf = null
params.reference = null
params.targets = null
params.threads = null
params.mem = null
params.HapCUT2 = null
params.interval = null
params.outdir = null
params.help = false

// TODO: make help string
// Help message
helpMessage = """
Parameters:
--outdir            Desired path/name of folder to store output in.
""".stripIndent()

// Show help when needed
if (params.help){
    log.info helpMessage
        exit 0
}

// Make sure necessary input parameters are assigned.
assert params.bam_paths != null, 'Input parameter "bam_paths" cannot be unasigned.'
assert params.vcf != null, 'Input parameter "vcf" cannot be unasigned.'
assert params.HapCUT2!= null, 'Input parameter "HapCUT2" cannot be unasigned.'
assert params.interval != null, 'Input parameter "interval" cannot be unasigned.'
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.targets != null, 'Input parameter "targets" cannot be unasigned.'
assert params.threads != null, 'Input parameter "threads" cannot be unasigned.'
assert params.mem != null, 'Input parameter "mem" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "bam_paths          : ${params.bam_paths}"
println "vcf                : ${params.vcf}"
println "HapCUT2            : ${params.HapCUT2}"
println "interval           : ${params.interval}"
println "reference          : ${params.reference}"
println "targets            : ${params.targets}"
println "threads            : ${params.threads}"
println "mem                : ${params.mem}"
println "outdir             : ${params.outdir}"

// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
bam_paths_ch = Channel.fromPath(params.bam_paths)
bam_paths_ch
    .splitCsv(header: true)
    .map { it -> [it.sample, file(it.bam_path), file(it.bai_path)] }
    .into { bam_paths_process_ch; bam_paths_check_ch }

// Map this channel to just the sample names.
bam_sample_names_ch = bam_paths_check_ch.map { it -> it[0] }

// Get file handlers for input files.
vcf = file(params.vcf)
reference = file(params.reference)  // Directory of 10x reference.
reference_fa = file(params.reference + '/fasta/genome.fa')  // Reference fasta file.
targets = file(params.targets)

// TODO: this is temporary, HapCUT2 should be installed properly, and this should not be necessary.
HAPCUT2 = params.HapCUT2 + "/build/HAPCUT2"
extractHAIRS = params.HapCUT2 + "/build/extractHAIRS"
LinkFragments = params.HapCUT2 + "/utilities/LinkFragments.py"

// FIXME:
// this is just for testing
process make_small_bam {
    input:
    set sample, file(bam), file(bai) from bam_paths_process_ch

    output:
    set sample, file("small.bam"), file("small.bam.bai") into small_bam_ch, small_bam_haplotag_bam_ch

    script:
    """
    samtools view -b -o "small.bam" -T $reference_fa $bam $params.interval
    samtools index -b "small.bam"
    """
}

// Extract sample names from VCF.
process get_sample_names {
    output:
    file "samples.txt" into vcf_sample_names_ch

    script:
    """
    # Get the line of the VCF containing the sample names.
    grep -m 1 "^#CHROM" $vcf > "header_line.txt"
    # Get the sample names in one line.
    cat "header_line.txt" | awk '{print substr(\$0, index(\$0,\$10))}' > "samples_line.txt"
    # Make a file with one line per sample name.
    cat "samples_line.txt" | tr -s '\t' '\n' > "samples.txt"
    """
}

// Convert this channel from a file with sample names, to a channel with one item per sample name.
vcf_sample_names_ch
    .splitText()
    .map { it.trim() }
    .into { vcf_sample_names_check_ch; vcf_sample_names_split_ch }

//vcf_sample_names_check_ch.subscribe { println it }
//vcf_sn = vcf_sample_names_check_ch.collect()
//bam_sn = bam_sample_names_ch.toList().toSet()
//println "${vcf_sn.getClass()}"
//print "${vcf_sn}"
//intersect_size = vcf_sn.intersect(bam_sn).size()
//union_size = vcf_sn.union(bam_sn).size()
//assert vcf_sn == bam_sn, "error"

// FIXME: -L option only for testing
// Extract a single sample from the VCF. This process is run once for each sample, splitting the
// VCF into as many files as there are samples.
// --remove-unused-alternates is used to avoid gentypes like "2/3", which HapCUT2 can't deal with.
// --exclude-non-variants is used to exclude missing genotypes, formatted as "./.", which HapCUT2
// can't deal with either. These two arguments don't take anything away from the data, as when
// the samples are merged together, everything will be the same again.
process get_sample_vcf {
    input:
    val sample from vcf_sample_names_split_ch

    output:
    //set sample, file("sample.vcf") into vcf_extract_ch, vcf_link_ch, vcf_phase_ch
    set sample, file("sample.vcf"), file("sample.vcf.idx") into vcf_ch

    script:
    """
    mkdir tmp
    gatk SelectVariants \
        -V $vcf \
        -R $reference_fa \
        -sn $sample \
        -L ${params.interval} \
        --remove-unused-alternates \
        -O "sample.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${params.mem}g -Xms${params.mem}g"
    """
}

// FIXME:
// this is temporary, another solution needed.
process remove_nocalls {
    input:
    set sample, file(vcf), file(idx) from vcf_ch

    output:
    set sample, file("nocalls_removed.vcf"), file("nocalls_removed.vcf.idx") into vcf_nocalls_removed_ch

    script:
    """
    # The grep pattern below searches for two different patterns: ./. and .|.
    grep -v "\\./\\.\\|\\.|\\." $vcf > "nocalls_removed.vcf"
    gatk IndexFeatureFile -F "nocalls_removed.vcf"
    """
}
vcf_ch = vcf_nocalls_removed_ch

// Make a channel with (sample ID, VCF, VCF index, BAM, BAM index) tuples.
vcf_ch.join(small_bam_ch).into { data_extract_ch; data_link_ch; data_phase_ch }

// Convert BAM file to the compact fragment file format containing only haplotype-relevant information.
process extract_hairs {
    input:
    set sample, file(vcf), file(idx), file(bam), file(bai) from data_extract_ch

    output:
    set sample, file("unlinked_fragment") into unlinked_fragments_ch

    script:
    """
    $extractHAIRS --10X 1 --region $params.interval --bam $bam --VCF $vcf --out "unlinked_fragment"
    """
}

// Add unlinked fragment file to tuple for next process.
data_link_ch = data_link_ch.join(unlinked_fragments_ch)

// Use LinkFragments to link fragments into barcoded molecules.
process link_fragments {
    input:
    set sample, file(vcf), file(idx), file(bam), file(bai), file(unlinked_fragments) from data_link_ch

    output:
    set sample, file("linked_fragments") into linked_fragments_ch

    script:
    """
    python3 $LinkFragments --bam $bam --VCF $vcf --fragments $unlinked_fragments --out "linked_fragments"
    """
}

// Add linked fragment file to tuple for next process.
data_phase_ch = data_phase_ch.join(linked_fragments_ch)

// Use HAPCUT2 to assemble fragment file into haplotype blocks.
process phase_vcf {
    input:
    set sample, file(vcf), file(idx), file(bam), file(bai), file(linked_fragments) from data_phase_ch

    output:
    file "haplotypes" into haplotypes_ch
    set sample, file("haplotypes.phased.VCF") into phased_vcf_ch

    script:
    """
    $HAPCUT2 --nf 1 --outvcf 1 --fragments $linked_fragments --VCF $vcf --output "haplotypes"
    """
}


// Compress and index VCF.
process index_and_zip_vcf {
    input:
    set sample, file(vcf) from phased_vcf_ch

    output:
    set sample, file("phased.vcf.gz"), file("phased.vcf.gz.tbi") into phased_vcf_haplotag_ch
    file "phased.vcf.gz" into phased_vcf_merge_ch

    script:
    """
    bgzip -c $vcf > "phased.vcf.gz"
    tabix "phased.vcf.gz"
    """
}

// Merge VCFs into a multi-sample VCF, and compress and index it.
process merge_phased_vcf {
    publishDir "${params.outdir}/phased_vcf", mode: 'copy'

    input:
    val vcf_list from phased_vcf_merge_ch.toList()

    output:
    //set file("phased_merged.vcf.gz"), file("phased_merged.vcf.gz.tbi") into phased_merged_ch
    set file("phased_merged.vcf.gz") into phased_merged_ch

    script:
    // Prepare input string for vcf-merge.
    vcf_list_str = (vcf_list as List)
        .join(' ')  // Join paths in single string.
    """
    #vcf-merge --ref-for-missing 0 $vcf_list_str | bgzip -c > "phased_merged.vcf.gz"
    vcf-merge --ref-for-missing 0 $vcf_list_str > "phased_merged.vcf"
    bgzip "phased_merged.vcf" > "phased_merged.vcf.gz"
    tabix "phased_merged.vcf.gz"
    """
}

// Take the compressed and indexed VCF from above and join by sample name with BAM files.
data_haplotag_ch = phased_vcf_haplotag_ch.join(small_bam_haplotag_bam_ch)

// Add haplotype information to BAM, tagging each read with a haplotype (when possible), using
// the haplotype information from the phased VCF.
process haplotag_bam {
    input:
    set sample, file(phased_vcf), file(idx), file(bam), file(bai) from data_haplotag_ch

    output:
    set sample, file("phased.bam") into phased_bam_ch

    script:
    """
    whatshap haplotag --ignore-read-groups --reference $reference_fa -o "phased.bam" $phased_vcf $bam
    """
}

process index_bam {
    publishDir "${params.outdir}/phased_bam/$sample", mode: 'copy'

    input:
    set sample, file(bam) from phased_bam_ch

    output:
    set sample, file("phased.bam"), file("phased.bam.bai") into indexed_phased_bam

    script:
    """
    samtools index "phased.bam"
    """
}
