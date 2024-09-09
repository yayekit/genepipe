import os
from pathlib import Path
import snakemake

def ngs_data_pipeline(input_dir, output_dir):
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Define Snakemake workflow
    snakefile = """
    import os

    INPUT_DIR = "{input_dir}"
    OUTPUT_DIR = "{output_dir}"
    REFERENCE_GENOME = "path/to/reference_genome.fa"

    # Get all input FASTQ files
    SAMPLES = [f.split('.')[0] for f in os.listdir(INPUT_DIR) if f.endswith('.fastq')]

    rule all:
        input:
            expand(OUTPUT_DIR + "{{sample}}_variants.vcf", sample=SAMPLES),
            OUTPUT_DIR + "multiqc_report.html",
            OUTPUT_DIR + "pipeline_metrics.txt"

    rule fastqc:
        input:
            INPUT_DIR + "{{sample}}.fastq"
        output:
            OUTPUT_DIR + "fastqc/{{sample}}_fastqc.html",
            OUTPUT_DIR + "fastqc/{{sample}}_fastqc.zip"
        shell:
            "fastqc {{input}} -o {OUTPUT_DIR}/fastqc"

    rule trimmomatic:
        input:
            INPUT_DIR + "{{sample}}.fastq"
        output:
            trimmed = OUTPUT_DIR + "trimmed/{{sample}}_trimmed.fastq",
            trimlog = OUTPUT_DIR + "trimmed/{{sample}}_trimlog.txt"
        shell:
            "trimmomatic SE {{input}} {{output.trimmed}} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 "
            "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> {{output.trimlog}}"

    rule post_trim_fastqc:
        input:
            OUTPUT_DIR + "trimmed/{{sample}}_trimmed.fastq"
        output:
            OUTPUT_DIR + "fastqc_post_trim/{{sample}}_trimmed_fastqc.html",
            OUTPUT_DIR + "fastqc_post_trim/{{sample}}_trimmed_fastqc.zip"
        shell:
            "fastqc {{input}} -o {OUTPUT_DIR}/fastqc_post_trim"

    rule bwa_mem:
        input:
            OUTPUT_DIR + "trimmed/{{sample}}_trimmed.fastq"
        output:
            OUTPUT_DIR + "aligned/{{sample}}_aligned.sam"
        shell:
            "bwa mem {REFERENCE_GENOME} {{input}} > {{output}}"

    rule sam_to_bam:
        input:
            OUTPUT_DIR + "aligned/{{sample}}_aligned.sam"
        output:
            OUTPUT_DIR + "aligned/{{sample}}.bam"
        shell:
            "samtools view -bS {{input}} > {{output}}"

    rule sort_bam:
        input:
            OUTPUT_DIR + "aligned/{{sample}}.bam"
        output:
            OUTPUT_DIR + "sorted/{{sample}}_sorted.bam"
        shell:
            "samtools sort {{input}} -o {{output}}"

    rule index_bam:
        input:
            OUTPUT_DIR + "sorted/{{sample}}_sorted.bam"
        output:
            OUTPUT_DIR + "sorted/{{sample}}_sorted.bam.bai"
        shell:
            "samtools index {{input}}"

    rule mark_duplicates:
        input:
            OUTPUT_DIR + "sorted/{{sample}}_sorted.bam"
        output:
            bam = OUTPUT_DIR + "dedup/{{sample}}_dedup.bam",
            metrics = OUTPUT_DIR + "dedup/{{sample}}_dedup_metrics.txt"
        shell:
            "picard MarkDuplicates I={{input}} O={{output.bam}} M={{output.metrics}}"

    rule base_recalibration:
        input:
            bam = OUTPUT_DIR + "dedup/{{sample}}_dedup.bam",
            known_sites = "path/to/known_sites.vcf"
        output:
            recal_table = OUTPUT_DIR + "recal/{{sample}}_recal_data.table"
        shell:
            "gatk BaseRecalibrator -I {{input.bam}} -R {REFERENCE_GENOME} "
            "--known-sites {{input.known_sites}} -O {{output.recal_table}}"

    rule apply_bqsr:
        input:
            bam = OUTPUT_DIR + "dedup/{{sample}}_dedup.bam",
            recal_table = OUTPUT_DIR + "recal/{{sample}}_recal_data.table"
        output:
            OUTPUT_DIR + "recal/{{sample}}_recal.bam"
        shell:
            "gatk ApplyBQSR -I {{input.bam}} -R {REFERENCE_GENOME} "
            "--bqsr-recal-file {{input.recal_table}} -O {{output}}"

    rule variant_calling:
        input:
            bam = OUTPUT_DIR + "recal/{{sample}}_recal.bam"
        output:
            OUTPUT_DIR + "variants/{{sample}}_variants.vcf"
        shell:
            "gatk HaplotypeCaller -R {REFERENCE_GENOME} -I {{input.bam}} -O {{output}}"

    rule annotate_variants:
        input:
            OUTPUT_DIR + "variants/{{sample}}_variants.vcf"
        output:
            OUTPUT_DIR + "variants/{{sample}}_annotated.vcf"
        shell:
            "snpEff -v GRCh38.86 {{input}} > {{output}}"

    rule filter_variants:
        input:
            OUTPUT_DIR + "variants/{{sample}}_annotated.vcf"
        output:
            OUTPUT_DIR + "variants/{{sample}}_filtered.vcf"
        shell:
            "bcftools filter -i 'QUAL>20 && DP>10' {input} | bcftools filter -e 'INFO/DP<10 || INFO/DP>100' > {output}"

    rule combine_variants:
        input:
            expand(OUTPUT_DIR + "variants/{{sample}}_filtered.vcf", sample=SAMPLES)
        output:
            OUTPUT_DIR + "variants/combined_filtered_variants.vcf"
        shell:
            "bcftools merge {{input}} -o {{output}}"

    rule variant_summary_plot:
        input:
            OUTPUT_DIR + "variants/combined_filtered_variants.vcf"
        output:
            OUTPUT_DIR + "plots/variant_summary.png"
        script:
            "scripts/plot_variant_summary.py"

    rule sample_similarity_plot:
        input:
            OUTPUT_DIR + "variants/combined_filtered_variants.vcf"
        output:
            OUTPUT_DIR + "plots/sample_similarity.png"
        script:
            "scripts/plot_sample_similarity.py"

    rule multiqc:
        input:
            expand(OUTPUT_DIR + "fastqc/{sample}_fastqc.zip", sample=SAMPLES),
            expand(OUTPUT_DIR + "fastqc_post_trim/{sample}_trimmed_fastqc.zip", sample=SAMPLES),
            expand(OUTPUT_DIR + "dedup/{sample}_dedup_metrics.txt", sample=SAMPLES),
            expand(OUTPUT_DIR + "variants/{sample}_variants.vcf", sample=SAMPLES)
        output:
            OUTPUT_DIR + "multiqc_report.html"
        shell:
            "multiqc {OUTPUT_DIR} -o {OUTPUT_DIR}"

    rule collect_metrics:
        input:
            expand(OUTPUT_DIR + "variants/{sample}_variants.vcf", sample=SAMPLES),
            OUTPUT_DIR + "multiqc_report.html"
        output:
            OUTPUT_DIR + "pipeline_metrics.txt"
        run:
            import glob
            import subprocess

            with open(output[0], "w") as metrics_file:
                # Count total reads
                total_reads = sum(int(subprocess.check_output(f"zcat {{INPUT_DIR}}/{sample}.fastq.gz | wc -l", shell=True)) // 4 for sample in SAMPLES)
                metrics_file.write(f"Total reads: {{total_reads}}\\n")

                # Count reads after trimming
                trimmed_reads = sum(int(subprocess.check_output(f"zcat {{OUTPUT_DIR}}/trimmed/{sample}_trimmed.fastq.gz | wc -l", shell=True)) // 4 for sample in SAMPLES)
                metrics_file.write(f"Reads after trimming: {{trimmed_reads}}\\n")

                # Calculate alignment rate
                for sample in SAMPLES:
                    aligned_reads = int(subprocess.check_output(f"samtools view -c {{OUTPUT_DIR}}/sorted/{{sample}}_sorted.bam", shell=True))
                    alignment_rate = (aligned_reads / trimmed_reads) * 100
                    metrics_file.write(f"Alignment rate for {{sample}}: {{alignment_rate:.2f}}%\\n")

                # Count variants
                for sample in SAMPLES:
                    variant_count = int(subprocess.check_output(f"grep -v '^#' {{OUTPUT_DIR}}/variants/{{sample}}_variants.vcf | wc -l", shell=True))
                    metrics_file.write(f"Variants called for {{sample}}: {{variant_count}}\\n")

    """

    # Write Snakefile
    with open("Snakefile", "w") as f:
        f.write(snakefile.format(input_dir=input_dir, output_dir=output_dir))

    # Run Snakemake
    snakemake.snakemake("Snakefile", cores=1, printshellcmds=True, dryrun=False)

    print("NGS data processing pipeline completed.")

# Usage example
# ngs_data_pipeline("/path/to/input/directory", "/path/to/output/directory")
