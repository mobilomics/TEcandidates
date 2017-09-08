![tecandidates_logos](https://user-images.githubusercontent.com/31257257/29850516-b42a2540-8d04-11e7-9286-7e8ebefbd455.png)

# TEcandidates
TEcandidates is a pipeline to include transposable elements in RNA-seq differential expression analysis.


### INSTALLATION INSTRUCTIONS
_________________________________________________

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
_________________________________________________

After installing the required softwares, grant execution permissions to the pipeline script:

    $ chmod u+x TEcandidates.sh

To run TEcandidates, the script must be executed as:

    $ ./TEcandidates.sh -t=Number_of_threads -r=RAM_to_use -g=Genome_Fasta_File -fq=Path_to_FASTQ_files -m=Mode -c=Coverage -te=TE_Annotation

    -t Number of threads to use in the softwares executed during the pipeline
    -r Maximum amount of RAM assigned to Trinity (Trinity's --max_memory option)
    -g Genome to use (FASTA format, .fasta extension)
    -fq Path to FASTQ files to use (all files must have .fastq extension)
    -m Mode of FASTQ files, SE for Single-end reads and PE for Paired-end reads
    -c Minimum coverage in which a Transposable element must be covered by a de-novo transcript in order to be selected as candidate


TEcandidates can be used with either single-end reads or paired-end reads. Reads files must have ".fastq" extension.

### SAMPLE DATA
_________________________________________________

In order to check that TEcandidates is working correctly, please test it with a dataset from _Drosophila melanogaster_ (Ohtani et al., 2013). The dataset is available at Gene Expression Omnibus, accession no. GSE47006, and must be downloaded with The _fastq-dump_ tool from the SRA toolkit . To install the SRA toolkit, please copy the link of the appropriate version for your Operating System from https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/, and download it:

    $ wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

Check if fastq-dump is working:

    $ tar -xvzf sratoolkit.current-ubuntu64.tar.gz
    $ ./sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump -V

      ./sratoolkit.2.8.2-1-ubuntu64/bin/fastq-dump : 2.8.2

For ease of use, add the SRA toolkit to your $PATH environment variable. Get the full path to the SRA toolkit binaries with 

    $ readlink -f sratoolkit.2.8.2-1-ubuntu64/bin/

and append it to $PATH like this:

    $ export PATH=$PATH:path_to_SRAtoolkit


Create a new work directory

    $ mkdir TEcandidates_test
    $ cd TEcandidates_test

Download the _Drosophila melanogaster_ control dataset with:

    $ nohup fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR851837 > SRR851837.fastq-dump.log &

Download the _Drosophila melanogaster_ treatment dataset with:

    $ nohup fastq-dump --defline-seq '@$sn[_$rn]/$ri' --split-files --accession SRR851838 > SRR851838.fastq-dump.log &
    
Once these processes are done, check that they were downloaded correctly with 

    $ tail *log
    ==> SRR851837.fastq-dump.log <==
    Read 42134407 spots for SRR851837
    Written 42134407 spots for SRR851837

    ==> SRR851838.fastq-dump.log <==
    Read 48277060 spots for SRR851838
    Written 48277060 spots for SRR851838


Execute the pipeline script afterwards:






### CONTACT
_________________________________________________


### REFERENCES
_________________________________________________

Ohtani H, Iwasaki YW, Shibuya A, Siomi H, Siomi MC, Saito K. (2013). DmGTSF1 is necessary for Piwi-piRISC-mediated transcriptional transposon silencing in the Drosophila ovary. Genes Dev. 2013 Aug 1;27(15):1656-61. doi: 10.1101/gad.221515.113.



