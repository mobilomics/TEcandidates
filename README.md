![tecandidates_logos](https://user-images.githubusercontent.com/31257257/29850516-b42a2540-8d04-11e7-9286-7e8ebefbd455.png)

# TEcandidates
TEcandidates is a pipeline to include transposable elements in RNA-seq differential expression analysis.


### INSTALLATION INSTRUCTIONS

TEcandidates is implemented in Bash, and requires no installation. However, other softwares that are part of the pipeline are required:

**-BEDtools v2.25**

Download from https://github.com/arq5x/bedtools2/releases/tag/v2.25.0

    $ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
    $ tar -xvzf bedtools-2.25.0.tar.gz
    $ cd bedtools2/
    $ make
    $ make install

Check if installation was correct using        

    $ bedtools --version
    bedtools v2.25.0

**-BioPerl**

Please follow the instructions at http://bioperl.org/INSTALL.html

        
**-Bowtie v2.3**

Download from https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2

    $ unzip bowtie2-2.3.2-linux-x86_64.zip
    $ cd bowtie2-2.3.2/
    $ pwd

The above will show the full path to the bowtie2 binaries. Copy it and add it to the $PATH environment variable:

    $ export PATH=$PATH:path_to_bowtie2-2.3.2

Check correct installation with

    $ bowtie2 --version
    bowtie2-align-s version 2.3.2
    64-bit
    Built on dde45b53bd81
    Sat May  6 02:39:52 UTC 2017
    Compiler: gcc version 4.1.2 20080704 (Red Hat 4.1.2-55)
    Options: -O3 -m64 -msse2 -funroll-loops -g3 -DPOPCNT_CAPABILITY -DWITH_TBB -DNO_SPINLOCK -DWITH_QUEUELOCK=1
    Sizeof {int, long, long long, void*, size_t, off_t}: {4, 8, 8, 8, 8, 8}


**-Samtools v1.4.1**

Please follow the instructions at https://github.com/samtools/samtools


**-Trinity v2.4**

Download from https://github.com/trinityrnaseq/trinityrnaseq/

    $ wget https://github.com/trinityrnaseq/trinityrnaseq/archive/master.zip
    $ unzip master.zip
    $ cd trinityrnaseq-master/
    $ make
    $ pwd 

Copy the path to the Trinity directory and add it to the $PATH environment variable:

    $ export PATH=$PATH:path_to_Trinity

Check correct installation with

    $ Trinity --version
    Trinity version: Trinity-v2.4.0

### SAMPLE USAGE

After installing the required softwares, grant execution permissions to the pipeline script:

    $ chmod u+x TEcandidates.sh

To run TEcandidates, the script must be executed as:

    $ ./TEcandidates.sh -t=Number_of_threads -r=RAM_to_use -g=Genome_Fasta_File -fq=Path_to_FASTQ_files -m=Mode -c=Coverage -te=TE_Annotation

    -t Number of threads to use in the softwares executed during the pipeline
    -r Maximum amount of RAM assigned to Trinity (Trinity's --max\_memory option)
    -g Genome to use (FASTA format, .fasta extension)
    -fq Path to FASTQ files to use (all files must have .fastq extension)
    -m Mode of FASTQ files, SE for Single-end reads and PE for Paired-end reads
    -c Minimum coverage in which a Transposable element must be covered by a de-novo transcript in order to be selected as candidate


TEcandidates can be used with either single-end reads or paired-end reads. Reads files must have ".fastq" extension.



### CONTACT


### REFERENCES
