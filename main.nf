#!/usr/bin/env nextflow
/*
Author: Ólavur Mortensen <olavur@fargen.fo>
*/

// Input parameters.
params.bam_paths = null
params.vcf = null
params.reference = null
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
assert params.reference != null, 'Input parameter "reference" cannot be unasigned.'
assert params.outdir != null, 'Input parameter "outdir" cannot be unasigned.'

println "P I P E L I N E     I P U T S    "
println "================================="
println "bam_paths          : ${params.bam_paths}"
println "vcf                : ${params.vcf}"
println "reference          : ${params.reference}"
println "outdir             : ${params.outdir}"

// Turn the file with FASTQ paths into a channel with [sample, path] tuples.
bam_paths_ch = Channel.fromPath(params.bam_paths)
bam_paths_ch
    .splitCsv(header: true)
    .map { it -> [it.sample, file(it.bam_path), file(it.bai_path)] }
    .into { bam_paths_process_ch; bam_paths_check_ch }

// Map this channel to just the sample names.
bam_sample_names_ch = bam_paths_check_ch.map { it -> it[0] }

// Input multi sample VCF file.
vcf = file(params.vcf)

// Reference fasta file.
reference = file(params.reference)

// Extract a single sample from the VCF. This process is run once for each sample, splitting the
// VCF into as many files as there are samples.
// --remove-unused-alternates is used to avoid gentypes like "2/3", which HapCUT2 can't deal with.
// --exclude-non-variants is used to exclude missing genotypes, formatted as "./.", which HapCUT2
// can't deal with either. These two arguments don't take anything away from the data, as when
// the samples are merged together, everything will be the same again.
process get_sample_vcf {
    input:
    val sample from bam_sample_names_ch

    output:
    //set sample, file("sample.vcf") into vcf_extract_ch, vcf_link_ch, vcf_phase_ch
    set sample, file("sample.vcf"), file("sample.vcf.idx") into vcf_sample_ch

    script:
    """
    mkdir tmp
    gatk SelectVariants \
        -V $vcf \
        -R $reference \
        -sn $sample \
        --remove-unused-alternates \
        -O "sample.vcf" \
        --tmp-dir=tmp \
        --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g"
    """
}

// Remove the "PP" field from the VCF.
// PP is the posterior probability of the possible genotypes, calculated by GATK's CalculateGenotypePosteriors.
// VCFtools doesn't know how to merge this field, and raises an error, therefore we remove it here.
process reformat_vcf_pp {
    input:
    set sample, file(vcf), file(idx) from vcf_sample_ch

    output:
    set sample, file("reformat_pp.vcf") into vcf_reformat_pp_ch

    script:
    """
    remove_format_field.py --vcf $vcf --field PP > reformat_pp.vcf
    """
}

// Remove the "PL" field, for a similar reason as the "PP" field above.
// NOTE:
// The VCF contains sites where there is a reference allele and a missing alternate allele ("."). This is because
// the genotype is homozygous reference. In these cases, the PL field contains only one likelihood, because there
// is only one possible genotype. However, VCFtools insists there are two alleles (even though one of them is
// missing), and so there are three possible genotypes and therefore three likelihoods.
process reformat_vcf_pl {
    input:
    set sample, file(vcf) from vcf_reformat_pp_ch

    output:
    set sample, file("reformat_pl.vcf") into vcf_reformat_pl_ch

    script:
    """
    remove_format_field.py --vcf $vcf --field PL > reformat_pl.vcf
    """
}

// Remove no-call variants (coded as "./.") from the VCF, as HapCUT2 doesn't seem to
// be able to handle these.
process remove_nocalls {
    input:
    set sample, file(vcf) from vcf_reformat_pl_ch

    output:
    set sample, file("nocalls_removed.vcf") into vcf_nocalls_removed_ch

    script:
    """
    remove_nocall_sites.py --vcf $vcf > nocalls_removed.vcf
    """
}

// Make a channel with (sample ID, VCF, VCF index, BAM, BAM index) tuples.
vcf_nocalls_removed_ch.join(bam_paths_process_ch).into { data_extract_ch; data_link_ch; data_phase_ch }

// NOTE: HapCUT2 does not accept compressed VCFs.

// Convert BAM file to the compact fragment file format containing only haplotype-relevant information.
process extract_hairs {
    input:
    set sample, file(vcf), file(bam), file(bai) from data_extract_ch

    output:
    set sample, file("unlinked_fragment") into unlinked_fragments_ch

    script:
    """
    extractHAIRS --10X 1 --bam $bam --VCF $vcf --out "unlinked_fragment"
    """
}

// Add unlinked fragment file to tuple for next process.
data_link_ch = data_link_ch.join(unlinked_fragments_ch)

// Use LinkFragments to link fragments into barcoded molecules.
process link_fragments {
    input:
    set sample, file(vcf), file(bam), file(bai), file(unlinked_fragments) from data_link_ch

    output:
    set sample, file("linked_fragments") into linked_fragments_ch

    script:
    """
    LinkFragments.py --bam $bam --VCF $vcf --fragments $unlinked_fragments --out "linked_fragments"
    """
}

// Add linked fragment file to tuple for next process.
data_phase_ch = data_phase_ch.join(linked_fragments_ch)

// Use HAPCUT2 to assemble fragment file into haplotype blocks.
// NOTE: this step may prune some low confidence genotypes, introducing no-calls into the dataset.
process phase_vcf {
    input:
    set sample, file(vcf), file(bam), file(bai), file(linked_fragments) from data_phase_ch

    output:
    file "haplotypes" into haplotypes_ch
    set sample, file("haplotypes.phased.VCF") into phased_vcf_ch

    script:
    """
    HAPCUT2 --nf 1 --outvcf 1 --fragments $linked_fragments --VCF $vcf --output "haplotypes"
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
// FIXME: Nextflow can't deal with "0\|0", because of the backslash, so I used a double backslash,
// but I'm not sure if this will work.
process merge_phased_vcf {
    publishDir "${params.outdir}/phased_vcf", mode: 'copy', overwrite: true

    input:
    val vcf_list from phased_vcf_merge_ch.toList()

    output:
    file("phased_merged.vcf.gz") into merged_vcf_ch

    script:
    // Prepare input string for vcf-merge.
    vcf_list_str = (vcf_list as List)
        .join(' ')  // Join paths in single string.
    """
    vcf-merge --ref-for-missing 0\\|0 $vcf_list_str > "phased_merged.vcf"
    bgzip -c "phased_merged.vcf" > "phased_merged.vcf.gz"
    """
}

// Index merged VCF.
process index_merged_vcf {
    publishDir "${params.outdir}/phased_vcf", mode: 'copy', pattern: '*.vcf.gz.tbi', overwrite: true

    input:
    file(vcf) from merged_vcf_ch

    output:
    file("*.vcf.gz.tbi")

    script:
    """
    tabix $vcf
    """
}
