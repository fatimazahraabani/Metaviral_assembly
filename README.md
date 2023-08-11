# Metaviral_assembly

## Description
This Snakemake workflow is designed to perform viral metagenomic assembly and analysis. It takes paired-end FASTQ files, performs mapping, assembly, optimization, and evaluation of the assembled contigs using various tools. The output includes optimized assemblies, QUAST reports, Diamond alignment results, ORF detection, RdRp detection, and chimera detection.

## Dependencies
- Snakemake
- Python >= 3.6
- BWA (0.7.17)
- Picard (2.23.5)
- Megahit (1.2.9)
- Spades (3.15.5)
- MetaSpades (3.15.5)
- MetaViralSpades (3.15.5)
- RNASpades (3.15.5)
- Cap3 (10.2011)
- QUAST (5.0.2)
- Diamond (latest version)
- HMMER (latest version)
- Samtools (installed as a dependency for other tools)

## Setup

1. Install Snakemake: Make sure you have Snakemake installed. If not, you can install it via pip:
   ```
   pip install snakemake
   ```

2. Install the required tools: Ensure that all the listed tools and their dependencies are installed and available in the system path.

3. Prepare input data: Place your paired-end FASTQ files in the specified `path_data` directory. The workflow expects paired-end reads named `{sample}_1.fastq` and `{sample}_2.fastq`.

4. Clone this repository:

```bash
git clone https://github.com/yourusername/viral-metagenomics-workflow.git
cd viral-metagenomics-workflow

## Configuration

Before running the workflow, you need to specify the configuration parameters in a JSON file (`cluster.json`). The configuration file should contain parameters for the number of threads, memory allocation, and other relevant cluster options required to run the tools efficiently on your computing infrastructure.

## Running the Workflow

To execute the workflow, navigate to the directory containing the Snakefile and the configuration file. Then, run the following command:

```
snakemake --use-conda
```

The `--use-conda` flag enables Snakemake to automatically create and manage the required Conda environments for the rules that specify the `envmodules` directive.

## Output

The workflow generates the following output:

- Assemblies for different assemblers (Megahit, Spades, MetaSpades, MetaViralSpades, RNASpades) in the results/ directory.
- QUAST evaluation reports for each assembly in the results/ directory.
- Diamond alignment results in the results/ directory.
- ORF detection results for each assembly.
- RdRp detection results for each assembly.
- RdRp HMM search results for each assembly.
- Chimera detection results for each assembly.

## Note

- The workflow assumes that the required reference genome (`GCF_016801865.2_TS_CPP_V2_genomic.fasta`) is available in the specified path. Make sure to provide the correct path to the reference genome in the `bwa_mapping` rule.
- The workflow expects specific naming conventions for the input FASTQ files and may require modifications if your files have different naming patterns.

## Contact Information

For any questions or issues related to this workflow, please contact abanifatimazahra@gmail.com

