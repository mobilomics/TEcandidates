![tecandidates_logos](https://user-images.githubusercontent.com/31257257/29850516-b42a2540-8d04-11e7-9286-7e8ebefbd455.png)

# TEcandidates
TEcandidates is a pipeline to include transposable elements in RNA-seq differential expression analysis.


### INSTALLATION INSTRUCTIONS

TEcandidates is implemented in Bash, and requires no installation. However, other softwares that are part of the pipeline are required:

**-BEDtools v2.25**

    def foo():
        if not bar:
            return True
        
**-Bowtie v2.3**


**-Perl v5.20.2**


**-Samtools v1.4.1**


**-Trinity v2.4**


### SAMPLE USAGE

After installing the required softwares, grant execution permissions to the pipeline script:

    chmod u+x TEcandidates.sh

To run TEcandidates, the script must be executed as:

    ./TEcandidates.sh -t=Number_of_threads -r=RAM_to_use -g=Genome_Fasta_File -fq=Path_to_FASTQ_files -m=Mode -c=Coverage -te=TE_Annotation

    -t Number of threads to use in the softwares executed during the pipeline
    -r Maximum amount of RAM assigned to Trinity (Trinity's --max\_memory option)
    -g Genome to use (FASTA format, .fasta extension)
    -fq Path to FASTQ files to use (all files must have .fastq extension)
    -m Mode of FASTQ files, SE for Single-end reads and PE for Paired-end reads
    -c Minimum coverage in which a Transposable element must be covered by a de-novo transcript in order to be selected as candidate


TEcandidates can be used with either single-end reads or paired-end reads. Reads files must have ".fastq" extension.
