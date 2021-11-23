/* 
 * pipeline input parameters 
 */
params.cwd = "/mnt/volume1/data"
params.reads_summary = "$baseDir/aws_data/ON001-RNA-R*.{fastq.gz,seq_summary_location}"
//params.reads = "$baseDir/*/*.fastq.gz"
//params.seq_summary = "$baseDir/*/*.seq_summary_location"
params.fast5_parent_folder = "$baseDir/aws_data"
params.multiqc = "$baseDir/multiqc"
params.outdir = "xpore_results"
params.ld_preload = "/home/ubuntu/miniconda3/envs/xpore/lib/libasan.so"

params.yml = "$baseDir/*.yml"

// from observation does not use more than 6 threads 
params.threads = 8

params.genome = "/mnt/volume1/genome/Macaca_fascicularis.Macaca_fascicularis_6.0.dna.primary_assembly.fa"
params.gtf = "/mnt/volume1/genome/Macaca_fascicularis.Macaca_fascicularis_6.0.104.gtf"
params.transcriptome = "/mnt/volume1/genome/transcriptome/Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa"
params.transcriptome_index = "/mnt/volume1/genome/transcriptome/mfa_6.0_cdna.mmi"

log.info """\
         R N A S E Q - N F   P I P E L I N E    
         ===================================
         Current Dir : ${params.cwd}
		 seq_summary : ${params.seq_summary}
         reads : ${params.reads}
		 fast5 : ${params.fast5_parent_folder}
         outdir : ${params.outdir}
         """
         .stripIndent()

 
/* 
* 	Create Channel
*	reads_ch = Channel
*			.fromFilePairs( params.reads )
*			.view()
*	reads_ch.into { reads_minimap2_ch ; reads_nanopolish_index_ch ; reads_nanopolish_eventalign_ch ; reads_xpore_dataprep_ch }
*	seqSum_ch 	= 	Channel
*				.fromFilePairs( params.seq_summary )
*/

reads_summary_ch = 	Channel
					.fromFilePairs( params.reads_summary )

reads_summary_ch.into { reads_minimap2_ch ; reads_nanopolish_index_ch ; reads_nanopolish_eventalign_ch ; reads_xpore_dataprep_ch }
			

/* 
* 	https://xpore.readthedocs.io/en/latest/preparation.html
* 	Start with step 2 as basecalled data is already generated
*	Step 0_2: Align to transcriptome
*/
process minimap2_fastq {
    	
	tag "$sample_id"

    input:
    tuple val(sample_id) , path(fastq) from reads_minimap2_ch
	path transcript_index from params.transcriptome_index
	val threads from params.threads
	
    output:
	tuple path("${sample_id}.bam") , path("${sample_id}.bam.bai")  into minimap2_sam_ch
    
	script:   
    """
    minimap2 -ax map-ont -uf -t $threads --secondary=no ${transcript_index} ${fastq[0]} > ${sample_id}.sam 2>> log_minimap2.txt
	samtools view -@ $threads -Sb ${sample_id}.sam | samtools sort -o ${sample_id}.bam - 2>> log_sam2bam.txt
	samtools index ${sample_id}.bam >> log_bamIndex.txt
	rm ${sample_id}.sam
	"""
}


/* 
*	Step 0_3a: nanopolish index
*/
process nanopolish_index {
    	
	tag "$sample_id"

    input:
    tuple val(sample_id) , path(fastq) from reads_nanopolish_index_ch
	//tuple val(runID) , path(seq_summary) from seqSum_ch
	val fast5 from params.fast5_parent_folder
	
    output:
	path "${sample_id}.fastq.gz.index" into reads_index_ch
	path "${sample_id}.fastq.gz.index.fai" into reads_fai_ch
	path "${sample_id}.fastq.gz.index.gzi" into reads_gzi_ch
	path "${sample_id}.fastq.gz.index.readdb" into reads_readdb_ch
    
	script:   
    """
	nanopolish index -f ${fastq[1]} -d $fast5/$sample_id ${fastq[0]}
	"""
}

/* 
*	Step 0_3b: nanopolish eventalign
*/
process nanopolish_eventalign {
    
	validExitStatus 0,1	
	
	tag "$sample_id"
	publishDir path: "${params.outdir}/${sample_id}/0_nanopolish" , mode: 'copy' , pattern: '*.txt'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    tuple val(sample_id) , path(fastq) from reads_nanopolish_eventalign_ch
	tuple path(bam) , path(bam_bai) from minimap2_sam_ch
	val genome from params.transcriptome
	val ld_preload from params.ld_preload
	path(index) from reads_index_ch
	path(index_fai) from reads_fai_ch
	path(index_gzi) from reads_gzi_ch
	path(index_readdb) from reads_readdb_ch
	val threads from params.threads
	
    output:
	file "${sample_id}_summary.txt" into summary_ch
	file "${sample_id}_eventalign.txt" into eventalign_ch
    file 'SUCCESS_*' into log_ch_0_3b
	val sample_id into sampleID_ch_0_3b
	
	script:   
    """
	LD_PRELOAD=$ld_preload \
	nanopolish eventalign \
	--reads ${fastq[0]} \
	--bam $bam \
	--genome $genome \
	--signal-index --scale-events --summary ${sample_id}_summary.txt \
	--threads 16 > ${sample_id}_eventalign.txt 2>> SUCCESS_step0_3b_Nanopolish_eventalign
	echo "COMPLETED Step0 (Nanopolish_eventalign) : ${sample_id}" >> SUCCESS_step0_3b_Nanopolish_eventalign
    """
}

/* 
*	Step 1: xpore dataprep
*/
process xpore_dataprep {
    	
	tag "$sample_id"
	publishDir path: "${params.outdir}/${sample_id}" , mode: 'copy'  , pattern: 'dataprep/*'
	publishDir path: "${params.outdir}/${sample_id}/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
	path eventalign from eventalign_ch
	path gtf from params.gtf
	val threads from params.threads
	val transcriptome from params.transcriptome
	val sample_id from sampleID_ch_0_3b
	
    output:
	file "dataprep/eventalign.index" into step1_eventalign_index_ch
	file "dataprep/data.json" into step1_data_json_ch
	file "dataprep/data.index" into step1_data_index_ch
    file "dataprep/data.readcount" into step1_data_readcount_ch
	file "dataprep/transcript_id_version_merged.gtf" into step1_trascript_merged_ch 
	file "dataprep/data.log" into step1_data_log_ch
	file 'SUCCESS_*' into log_ch_1
	val sample_id into sampleID_ch_1
	
	script:   
    """
	xpore dataprep --eventalign $eventalign \
	--gtf_path_or_url $gtf \
	--transcript_fasta_paths_or_urls $transcriptome \
	--genome \
	--out_dir dataprep \
	--n_processes $threads 2>> SUCCESS_step1_xpore_dataprep
	echo "COMPLETED Step1 (xPORE_dataprep) : ${sample_id}" >> SUCCESS_step1_xpore_dataprep
    """
}

/* 
*	Step 2: Maunal create your own YML file
*/

/* 
*	Step 3: xpore diffmod
*/
process xpore_diffmod {
    	
	tag "$sample_id"
	publishDir path: "${params.outdir}/_xpore-diffmod" , mode: 'copy'  , pattern: 'dataprep/*'
	publishDir path: "${params.outdir}/_xpore-diffmod/logs" , mode: 'copy' , pattern: 'SUCCESS_*'
	
    input:
    path comparison_yml from params.yml
	val sample_id from sampleID_ch_1.collect()
	
    output:
	file 'SUCCESS_*' into log_ch_3
	
	script:   
    """
	xpore diffmod --config ${comparison_yml} 2>> SUCCESS_step3_xpore_diffmod
	echo "COMPLETED Step3 (xPORE_diffmod) : ${sample_id}" >> SUCCESS_step3_xpore_diffmod
    """
}