# NGS Data Processing Pipeline

This project implements a Next-Generation Sequencing (NGS) data processing pipeline using Snakemake, a workflow management system. The pipeline automates the process of quality control, trimming, alignment, and post-alignment processing of NGS data.

## Features

- Quality control using FastQC
- Read trimming with Trimmomatic
- Alignment to a reference genome using BWA-MEM
- SAM to BAM conversion
- BAM file sorting and indexing

## Prerequisites

Before running this pipeline, ensure you have the following software installed:

- Python 3.6+
- Snakemake
- FastQC
- Trimmomatic
- BWA
- SAMtools

You can install the Python dependencies using pip:

```
pip install snakemake
```

For the other tools, please refer to their respective installation guides.

## Usage

1. Clone this repository:
   ```
   git clone https://github.com/yourusername/ngs-data-pipeline.git
   cd ngs-data-pipeline
   ```

2. Modify the `pipeline.py` script to set the correct paths for your input data, output directory, and reference genome.

3. Run the pipeline:
   ```
   python pipeline.py
   ```

   This will create a Snakefile and execute the workflow.

## Pipeline Steps

1. **FastQC**: Quality control check on raw sequence data.
2. **Trimmomatic**: Trim adapters and low-quality bases from reads.
3. **BWA-MEM**: Align trimmed reads to the reference genome.
4. **SAMtools**: Convert SAM to BAM, sort BAM files, and create BAM indexes.

## Customization

You can customize the pipeline by modifying the Snakemake rules in the `pipeline.py` file. Adjust parameters for each tool or add/remove steps as needed for your specific analysis.

## Output

The pipeline generates the following output files for each sample:

- `{sample}_fastqc.html`: FastQC report
- `{sample}_fastqc.zip`: Zipped FastQC results
- `{sample}_trimmed.fastq`: Trimmed reads
- `{sample}_aligned.sam`: Aligned reads in SAM format
- `{sample}.bam`: Aligned reads in BAM format
- `{sample}_sorted.bam`: Sorted BAM file
- `{sample}_sorted.bam.bai`: BAM index file

## Contributing

Contributions to improve the pipeline are welcome. Please fork the repository and submit a pull request with your changes.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- This pipeline uses several open-source bioinformatics tools. Please cite them appropriately in your research.
- Snakemake: Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2018.
