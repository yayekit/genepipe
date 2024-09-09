import os
import subprocess
from pathlib import Path

def ngs_data_pipeline(input_dir, output_dir):
    # Create output directory if it doesn't exist
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Step 1: Quality Control
    def run_fastqc(input_file, output_dir):
        cmd = f"fastqc {input_file} -o {output_dir}"
        subprocess.run(cmd, shell=True, check=True)

    # Step 2: Trimming
    def run_trimmomatic(input_file, output_dir):
        base_name = os.path.basename(input_file).split('.')[0]
        output_file = f"{output_dir}/{base_name}_trimmed.fastq"
        cmd = f"trimmomatic SE {input_file} {output_file} ILLUMINACLIP:TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
        subprocess.run(cmd, shell=True, check=True)
        return output_file

    # Step 3: Alignment
    def run_bwa(input_file, reference_genome, output_dir):
        base_name = os.path.basename(input_file).split('.')[0]
        output_file = f"{output_dir}/{base_name}_aligned.sam"
        cmd = f"bwa mem {reference_genome} {input_file} > {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        return output_file

    # Step 4: SAM to BAM conversion
    def sam_to_bam(input_file, output_dir):
        base_name = os.path.basename(input_file).split('.')[0]
        output_file = f"{output_dir}/{base_name}.bam"
        cmd = f"samtools view -bS {input_file} > {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        return output_file

    # Step 5: Sort BAM file
    def sort_bam(input_file, output_dir):
        base_name = os.path.basename(input_file).split('.')[0]
        output_file = f"{output_dir}/{base_name}_sorted.bam"
        cmd = f"samtools sort {input_file} -o {output_file}"
        subprocess.run(cmd, shell=True, check=True)
        return output_file

    # Step 6: Index BAM file
    def index_bam(input_file):
        cmd = f"samtools index {input_file}"
        subprocess.run(cmd, shell=True, check=True)

    # Main pipeline execution
    for file in os.listdir(input_dir):
        if file.endswith('.fastq'):
            input_file = os.path.join(input_dir, file)
            
            # Run FastQC
            run_fastqc(input_file, output_dir)
            
            # Run Trimmomatic
            trimmed_file = run_trimmomatic(input_file, output_dir)
            
            # Run BWA (assuming reference genome is available)
            reference_genome = "path/to/reference_genome.fa"
            aligned_file = run_bwa(trimmed_file, reference_genome, output_dir)
            
            # Convert SAM to BAM
            bam_file = sam_to_bam(aligned_file, output_dir)
            
            # Sort BAM file
            sorted_bam = sort_bam(bam_file, output_dir)
            
            # Index BAM file
            index_bam(sorted_bam)

    print("NGS data processing pipeline completed.")

# Usage example
# ngs_data_pipeline("/path/to/input/directory", "/path/to/output/directory")
