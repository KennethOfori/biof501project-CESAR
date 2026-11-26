// Resistant isolates: sample name -> SRA run accession. DOWNLOAD_READS fetches
// each run and names its FASTQs after the sample (e.g. PG1BP_1.fastq.gz).
params.samples = [
    PG1BP: 'SRR33692911',
    DMIH : 'SRR33692910',
    HCDL : 'SRR33692909',
    NTIH : 'SRR33692908'
]
params.reads_dir = "${projectDir}/data/reads"
params.outdir = "${projectDir}/results"
// K-12 reference as downloaded from NCBI Datasets; UNZIP_K12 extracts the files below.
// reference and gff are paths inside the archive, not on disk.
params.k12_zip   = "${projectDir}/data/k12.zip"
params.reference = "ncbi_dataset/data/GCF_000005845.2/GCF_000005845.2_ASM584v2_genomic.fna"
params.gff       = "ncbi_dataset/data/GCF_000005845.2/genomic.gff"

// Variant filtering. Bacteria are haploid, so a real variant should be carried by
// nearly every read at the site: min_af is deliberately strict.
params.min_qual = 30    // QUAL
params.min_dp   = 10    // read depth at the site
params.min_mq   = 30    // mapping quality, at pileup and as an INFO filter
params.min_bq   = 20    // base quality at pileup
params.min_af   = 0.9   // fraction of reads supporting the alt allele

// Extract the K-12 reference genome and annotation from the NCBI Datasets zip
process UNZIP_K12 {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'   // has bash + python; no unzip needed

    input:
    path zip

    output:
    path "k12/${params.reference}", emit: reference
    path "k12/${params.gff}",       emit: gff

    script:
    """
    python3 -m zipfile -e ${zip} k12
    """
}

// Download one SRA run and write its paired FASTQs under the sample name.
// storeDir keeps them in data/reads: if both files are already there, the
// download is skipped, even without -resume.
process DOWNLOAD_READS {
    label 'process_low'
    label 'error_retry'     // network hiccups
    tag "${sample_id} (${accession})"
    container 'quay.io/biocontainers/sra-tools:3.2.1--h4304569_1'
    storeDir params.reads_dir
    maxForks 2              // go easy on NCBI

    input:
    tuple val(sample_id), val(accession)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq.gz"), path("${sample_id}_2.fastq.gz")

    script:
    """
    prefetch ${accession} -O .
    fasterq-dump ./${accession} --split-files --threads ${task.cpus} --temp . -O .

    gzip -c ${accession}_1.fastq > ${sample_id}_1.fastq.gz
    gzip -c ${accession}_2.fastq > ${sample_id}_2.fastq.gz
    rm -rf ${accession} ${accession}_*.fastq
    """
}

// provide a quick QC summary of the fastq files
process FASTQC {
    label 'process_low'
    container 'quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0'
    publishDir "${params.outdir}/fastqc/", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc --threads ${task.cpus} -o . ${reads}

    """
}

// filter trim adapter sequences
process TRIMMOMATIC {
    label 'process_mid'
    tag "${sample_id}"
    container 'quay.io/biocontainers/trimmomatic:0.39--hdfd78af_2'
    publishDir "${params.outdir}/trimmomatic", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}_{1,2}.trimmed.fq.gz"), emit: reads
    path "${sample_id}_*.unpaired.fq.gz", emit: unpaired
    path "${sample_id}.trimmomatic.log", emit: log

    script:
    """
    # adapter FASTAs ship inside the bioconda package, next to the trimmomatic wrapper
    ADAPTERS=\$(dirname \$(readlink -f \$(which trimmomatic)))/adapters

    trimmomatic PE -threads ${task.cpus} \\
        ${reads[0]} ${reads[1]} \\
        ${sample_id}_1.trimmed.fq.gz ${sample_id}_1.unpaired.fq.gz \\
        ${sample_id}_2.trimmed.fq.gz ${sample_id}_2.unpaired.fq.gz \\
        ILLUMINACLIP:\${ADAPTERS}/TruSeq3-PE.fa:2:30:10:2:true \\
        LEADING:3 TRAILING:3 MINLEN:36 \\
        2> ${sample_id}.trimmomatic.log
    """
}

// Buiild the index file for BOWTIE2
process BOWTIE2_BUILD {
    label 'process_mid'
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'

    input:
    path fasta

    output:
    path "k12_index.*.bt2"

    script:
    """
    bowtie2-build --threads ${task.cpus} ${fasta} k12_index
    """
}
// align the reads to the K-12 E.coli strain using bowtie2
process ALIGNMENT_TO_REF {
    label 'process_mid'
    tag "${sample_id}"
    container 'quay.io/biocontainers/mulled-v2-ac74a7f02cebcfcc07d8e8d1d750af9c83b4d45a:f70b31a2db15c023d641c32f433fb02cd04df5a6-0'
    publishDir "${params.outdir}/k12_ref_alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}.sorted.bam"), path("${sample_id}.sorted.bam.bai"), emit: bam
    path "${sample_id}.bowtie2.log", emit: log

    script:
    """
    set -o pipefail
    bowtie2 -p ${task.cpus} -x k12_index -1 ${reads[0]} -2 ${reads[1]} 2> ${sample_id}.bowtie2.log |
        samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam -
    samtools index ${sample_id}.sorted.bam
    """
}

// report sample coverage after alignment
process SAMPLE_COVERAGE {
    label 'process_low'
    tag "${sample_id}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2'

    input:
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("${sample_id}.regions.bed.gz")

    script:
    """
    mosdepth -t ${task.cpus} -n --fast-mode --by 5000 ${sample_id} ${bam}
    """
}

// bcftools for variant calling (light weight processing)
process CALL_VARIANTS {
    label 'process_low'
    tag "${sample_id}"
    container 'quay.io/biocontainers/bcftools:1.25--h118bc1c_1'
    publishDir "${params.outdir}/variants", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam), path(bai)
    path reference

    output:
    tuple val(sample_id), path("${sample_id}.vcf.gz"), emit: vcf
    path "${sample_id}.raw.vcf.gz", emit: raw
    path "${sample_id}.filter_stats.tsv", emit: stats

    script:
    """
    set -o pipefail
    # -q/-Q discard poorly mapped reads and low-quality bases before calling
    bcftools mpileup -f ${reference} -q ${params.min_mq} -Q ${params.min_bq} -a FORMAT/AD ${bam} |
        bcftools call --ploidy 1 -mv -Oz -o ${sample_id}.raw.vcf.gz

    # Soft-filter first so the FILTER column records why a call was dropped, then
    # keep only PASS. -g/-G drop calls clustered around an indel, where the
    # alignment itself is unreliable.
    bcftools filter -g 10 -G 10 \\
        -e 'QUAL<${params.min_qual} || INFO/DP<${params.min_dp} || INFO/MQ<${params.min_mq} || (FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1])<${params.min_af}' \\
        -s LOWQUAL -Oz -o ${sample_id}.soft.vcf.gz ${sample_id}.raw.vcf.gz
    bcftools view -f PASS -Oz -o ${sample_id}.vcf.gz ${sample_id}.soft.vcf.gz

    printf 'sample\\traw_calls\\tpassed\\n' > ${sample_id}.filter_stats.tsv
    printf '%s\\t%s\\t%s\\n' ${sample_id} \\
        \$(bcftools view -H ${sample_id}.raw.vcf.gz | wc -l) \\
        \$(bcftools view -H ${sample_id}.vcf.gz | wc -l) >> ${sample_id}.filter_stats.tsv
    """

}

process CIRCOS_PLOT {
    label 'process_low'
    tag "${sample_id}"
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/circos", mode: 'copy'

    input:
    tuple val(sample_id), path(coverage), path(vcf)

    output:
    path "${sample_id}.circos.png"

    script:
    """
    circos_plot.py ${sample_id} ${coverage} ${vcf}
    """
}

// ---------------- Reference-based locus analysis ----------------
// Loci present in BOTH K-12 and the isolates, so a comparison is meaningful:
// porins, efflux regulators, drug targets. Explains phenotypes the acquired
// gene content does not, e.g. carbapenem resistance with no carbapenemase.

process BUILD_LOCUS_BED {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/loci", mode: 'copy'

    input:
    path gff

    output:
    path "resistance_loci.bed"

    script:
    """
    make_locus_bed.py ${gff} resistance_loci.bed
    """
}

// Per-locus depth; --by a BED instead of the 5 kb windows used for the circos
process LOCUS_COVERAGE {
    label 'process_low'
    tag "${sample_id}"
    container 'quay.io/biocontainers/mosdepth:0.3.14--h87be163_2'

    input:
    tuple val(sample_id), path(bam), path(bai)
    path bed

    output:
    tuple val(sample_id), path("${sample_id}.regions.bed.gz"), path("${sample_id}.mosdepth.summary.txt")

    script:
    """
    mosdepth -t ${task.cpus} -n --by ${bed} ${sample_id} ${bam}
    """
}

process LOCUS_TABLE {
    label 'process_low'
    tag "${sample_id}"
    container 'staphb/pycirclize:1.10.1'

    input:
    tuple val(sample_id), path(regions), path(summary), path(vcf)
    path bed

    output:
    path "${sample_id}.loci.tsv"

    script:
    """
    locus_table.py ${sample_id} ${bed} ${regions} ${summary} ${vcf} ${sample_id}.loci.tsv
    """
}

process LOCUS_HEATMAP {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/loci", mode: 'copy'

    input:
    path summary

    output:
    path "resistance_loci_heatmap.png"

    script:
    """
    locus_heatmap.py ${summary} resistance_loci_heatmap.png
    """
}

// ---------------- De novo assembly branch ----------------
// Recovers accessory genome (plasmids, mobile elements) that has no home in K-12

// create assembly 
process SPADES {
    label 'process_mid'
    tag "${sample_id}"
    container 'quay.io/biocontainers/spades:4.3.0--hde4eca7_2'
    publishDir "${params.outdir}/assembly", mode: 'copy', pattern: "*.{fasta,log}"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.fasta"), emit: contigs
    path "${sample_id}.spades.log", emit: log

    script:
    """
    spades.py --isolate \\
        -1 ${reads[0]} -2 ${reads[1]} \\
        -t ${task.cpus} -m ${task.memory.toGiga()} \\
        -o spades_out

    cp spades_out/contigs.fasta ${sample_id}.contigs.fasta
    cp spades_out/spades.log ${sample_id}.spades.log
    """
}

// Assembly QC for all samples in one report; K-12 given only for genome-fraction stats
process QUAST {
    label 'process_mid'
    container 'quay.io/biocontainers/quast:5.3.0--py313pl5321h5ca1c30_2'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path contigs
    path reference

    output:
    path "quast"

    script:
    """
    quast.py --threads ${task.cpus} -r ${reference} -o quast ${contigs}
    """
}

// Genome-wide relatedness between the isolates and against K-12. ANI settles
// "same species / same lineage?"; it saturates within a species, so it is not a
// substitute for core-genome SNP distances at finer scales.
process FAST_ANI {
    label 'process_mid'
    container 'quay.io/biocontainers/fastani:1.34--h4dfc31f_2'
    publishDir "${params.outdir}/ani", mode: 'copy'

    input:
    path contigs
    path reference

    output:
    path "ani.tsv"

    script:
    """
    ls ${contigs} > genomes.txt
    echo ${reference} >> genomes.txt
    fastANI --ql genomes.txt --rl genomes.txt -o ani.tsv -t ${task.cpus}
    """
}

process ANI_HEATMAP {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/ani", mode: 'copy'

    input:
    path ani

    output:
    path "ani_heatmap.png"

    script:
    """
    ani_heatmap.py ${ani} ani_heatmap.png
    """
}

// Acquired AMR genes + E. coli point mutations (gyrA, parC, ...). Image bundles the database.
process AMRFINDER {
    label 'process_low'
    tag "${sample_id}"
    container 'staphb/ncbi-amrfinderplus:4.2.7-2026-10-29.1'
    publishDir "${params.outdir}/amr", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)

    output:
    tuple val(sample_id), path("${sample_id}.amrfinder.tsv")

    script:
    """
    amrfinder --nucleotide ${contigs} \\
        --organism Escherichia --plus \\
        --name ${sample_id} \\
        --threads ${task.cpus} \\
        --output ${sample_id}.amrfinder.tsv
    """
}

// ---------------- Plasmid reconstruction ----------------
// MOB-suite's database is not in the image. Download it once; storeDir keeps it
// between runs so later runs skip this step.
process MOB_INIT {
    label 'process_low'
    container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1'
    storeDir "${projectDir}/data/mob_suite_db"

    output:
    path "mob_db"

    script:
    """
    mob_init -d mob_db
    """
}

// Label each contig as chromosome or plasmid and group plasmid contigs into plasmids
process MOB_RECON {
    label 'process_mid'
    tag "${sample_id}"
    container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1'
    publishDir "${params.outdir}/plasmids", mode: 'copy'

    input:
    tuple val(sample_id), path(contigs)
    path db

    output:
    // --prefix puts the sample id in front of every file MOB-recon writes
    tuple val(sample_id), path("${sample_id}/${sample_id}.contig_report.txt"), path("${sample_id}/${sample_id}.mobtyper_results.txt"), emit: reports
    path "${sample_id}/${sample_id}.*.fasta", emit: fasta

    script:
    """
    mob_recon -i ${contigs} -o ${sample_id} -d ${db} -n ${task.cpus} --prefix ${sample_id}

    # MOB-recon writes no mobtyper file when it finds no plasmids; keep the shape fixed
    [ -f ${sample_id}/${sample_id}.mobtyper_results.txt ] || touch ${sample_id}/${sample_id}.mobtyper_results.txt
    """
}

// Add chromosome/plasmid location to each AMRFinder hit
process AMR_LOCATION {
    label 'process_low'
    tag "${sample_id}"
    container 'quay.io/biocontainers/mob_suite:3.1.9--pyhdfd78af_1'
    publishDir "${params.outdir}/amr", mode: 'copy'

    input:
    tuple val(sample_id), path(amr), path(contig_report), path(mobtyper)

    output:
    path "${sample_id}.amr_location.tsv"

    script:
    """
    amr_plasmid_location.py ${sample_id} ${amr} ${contig_report} ${mobtyper} ${sample_id}.amr_location.tsv
    """
}

// Gene x sample grid, coloured by chromosome vs plasmid, for all samples together
process AMR_LOCATION_PLOT {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/amr", mode: 'copy'

    input:
    path summary

    output:
    path "amr_location_heatmap.png"

    script:
    """
    amr_location_heatmap.py ${summary} amr_location_heatmap.png
    """
}

// Linear maps of contigs carrying clustered resistance genes (integron / IS regions)
process RESISTANCE_REGION_PLOT {
    label 'process_low'
    container 'staphb/pycirclize:1.10.1'
    publishDir "${params.outdir}/amr", mode: 'copy'

    input:
    path summary

    output:
    path "resistance_regions.png"

    script:
    """
    resistance_region_plot.py ${summary} resistance_regions.png
    """
}

workflow {
    // emit a tupple pair [sample_id, [ read_1, read_2] ]
    accessions_ch = channel.fromList(params.samples.collect { name, acc -> tuple(name, acc) })
    reads_ch = DOWNLOAD_READS(accessions_ch).map { sample_id, r1, r2 -> tuple(sample_id, [r1, r2]) }

    // unzip the K-12 reference; both outputs are value channels, reusable by every sample
    UNZIP_K12(file(params.k12_zip, checkIfExists: true))
    ref = UNZIP_K12.out.reference
    gff = UNZIP_K12.out.gff

    FASTQC(reads_ch)
    TRIMMOMATIC(reads_ch)
    index_ch = BOWTIE2_BUILD(ref)
    ALIGNMENT_TO_REF(TRIMMOMATIC.out.reads, index_ch)

    SAMPLE_COVERAGE(ALIGNMENT_TO_REF.out.bam)
    CALL_VARIANTS(ALIGNMENT_TO_REF.out.bam, ref)
    CIRCOS_PLOT(SAMPLE_COVERAGE.out.join(CALL_VARIANTS.out.vcf))
    CALL_VARIANTS.out.stats.collectFile(
        name: 'variant_filter_stats.tsv',
        keepHeader: true,
        skip: 1,
        storeDir: "${params.outdir}/variants"
    )

    // per-locus coverage and variant burden against K-12
    loci_bed = BUILD_LOCUS_BED(gff)
    LOCUS_COVERAGE(ALIGNMENT_TO_REF.out.bam, loci_bed)
    LOCUS_TABLE(LOCUS_COVERAGE.out.join(CALL_VARIANTS.out.vcf), loci_bed)
    locus_summary = LOCUS_TABLE.out.collectFile(
        name: 'locus_summary.tsv',
        keepHeader: true,
        skip: 1,
        storeDir: "${params.outdir}/loci"
    )
    LOCUS_HEATMAP(locus_summary)

    // de novo assembly branch
    SPADES(TRIMMOMATIC.out.reads)
    QUAST(SPADES.out.contigs.map { _sample_id, contigs -> contigs }.collect(), ref)
    ANI_HEATMAP(FAST_ANI(SPADES.out.contigs.map { _sample_id, contigs -> contigs }.collect(), ref))
    AMRFINDER(SPADES.out.contigs)
    AMRFINDER.out.map { _sample_id, tsv -> tsv }.collectFile(
        name: 'amrfinder_summary.tsv',
        keepHeader: true,
        skip: 1,
        storeDir: "${params.outdir}/amr"
    )

    // plasmid reconstruction + AMR gene location
    MOB_RECON(SPADES.out.contigs, MOB_INIT())
    AMR_LOCATION(AMRFINDER.out.join(MOB_RECON.out.reports))
    location_summary = AMR_LOCATION.out.collectFile(
        name: 'amr_location_summary.tsv',
        keepHeader: true,
        skip: 1,
        storeDir: "${params.outdir}/amr"
    )
    AMR_LOCATION_PLOT(location_summary)
    RESISTANCE_REGION_PLOT(location_summary)
}



