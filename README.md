![tecandidates_logos](https://user-images.githubusercontent.com/31257257/29850516-b42a2540-8d04-11e7-9286-7e8ebefbd455.png)

# TEcandidates
TEcandidates is a pipeline to include transposable elements in RNA-seq differential expression analysis.


### INSTALLATION INSTRUCTIONS
_________________________________________________

#### 1. Dependencies

TEcandidates is implemented in Bash, and requires no installation. However, other softwares that are part of the pipeline are required. The following are the required softwares, along with some minimum instructions to install them in a computer with Linux. For additional help and/or troubleshooting, the source page of each software is listed for more detailed instructions.

- **BEDtools v2.25**

-Source page: https://github.com/arq5x/bedtools2/releases/tag/v2.25.0

    $ wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz
    $ tar -xvzf bedtools-2.25.0.tar.gz
    $ cd bedtools2/
    $ make
    $ sudo make install

Check if installation was correct using        

    $ bedtools --version
    bedtools v2.25.0

- **BioPerl**

-Source page: http://bioperl.org/INSTALL.html

		$ sudo perl -e shell -MCPAN
		Password:
		Terminal does not support AddHistory.

		cpan shell -- CPAN exploration and modules installation (v1.9800)
		Enter 'h' for help.

		cpan[1]> install C/CJ/CJFIELDS/BioPerl-1.007001.tar.gz

Verify installation with:

		$ perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"'

        
- **Bowtie v2.3**

-Source page: https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2

Download the appropiate file for your distribution and follow these instructions:

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


- **Samtools v1.4.1**

-Source page: https://github.com/samtools/samtools

		$ wget https://github.com/samtools/samtools/releases/download/1.4.1/samtools-1.4.1.tar.bz2
		$ tar -xvjf samtools-1.4.1.tar.bz2
		$ cd samtools-1.4.1
		$ ./configure
		$ make
		$ sudo make install

Check correct installation:
	
		$ samtools --version
		samtools 1.4.1
		Using htslib 1.4.1
		Copyright (C) 2017 Genome Research Ltd.

- **Trinity v2.4**

-Source page: https://github.com/trinityrnaseq/trinityrnaseq/

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

#### 2. Installing the TEcandidates pipeline

Download the TEcandidates tarball, and uncompress it:

    $ tar -xvzf TEcandidates_v1.tar.gz

Grant execution permissions to the pipeline script:

    $ chmod u+x TEcandidates_v1/TEcandidates.sh

For simplicity of use, add the TEcandidates full path to your PATH environment variable. First get the full path:

    $ readlink -f TEcandidates_v1
    
Then copy the output of the previous command, and add it to the PATH variable:

    $ export PATH=$PATH:/path/to/TEcandidates_v1

### SAMPLE USAGE
_________________________________________________

Once TEcandidates is in your PATH variable, you can execute it as

    $ TEcandidates.sh -t=Number_of_threads -r=RAM_to_use -g=Genome_Fasta_File -fq=Path_to_FASTQ_files -m=Mode -c=Coverage -te=TE_Annotation

    -t Number of threads to use in the softwares executed during the pipeline
    -r Maximum amount of RAM assigned to Trinity (Trinity's --max_memory option)
    -g Genome to use (FASTA format, .fasta extension)
    -fq Path to FASTQ files to use (all files must have .fastq extension)
    -m Mode of FASTQ files, SE for Single-end reads and PE for Paired-end reads
    -c Minimum coverage in which a Transposable element must be covered by a de-novo transcript in order to be selected as candidate

**Important considerations**
- Reads files must have ".fastq" extension.
- TEcandidates can be used with either single-end reads or paired-end reads.  
Paired-end reads **must have** "\_1.fastq" and "\_2.fastq" extensions.

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

    nohup TEcandidates.sh --threads=64 -g=dm3.fasta -fq=. -c=0.5 -te=dm3_rmsk_TE.gff3 -m=SE > TEcandidates.log &

Once it's done, you should have the following files:

    $ ls -lht
    total 28G
    drwxr-xr-x 3 user user 4.0K Sep 13 10:42 candidateTE_analysis_coverage-0.5
    -rw-r--r-- 1 user user 3.0M Sep 13 00:04 TEcandidates.log
    -rw------- 1 user user  328 Sep  6 12:43 nohup.out 
    -rw-r--r-- 1 user user 5.2M Sep  6 12:40 dm3_rmsk_TE.gff3
    -rw-r--r-- 1 user user 165M Sep  6 12:39 dm3.fasta 
    -rw-r--r-- 1 user user  15G Sep  6 07:26 SRR851838_1.fastq
    -rw-r--r-- 1 user user 3.8K Sep  6 07:26 SRR851838.fastq-dump.log
    -rw-r--r-- 1 user user  13G Sep  6 07:09 SRR851837_1.fastq
    -rw-r--r-- 1 user user 3.1K Sep  6 07:09 SRR851837.fastq-dump.log

The **candidateTE_analysis_coverage-0.5** folder contains the following files:

    $ ls -lht candidateTE_analysis_coverage-0.5/
    total 346M
    drwxr-xr-x 2 user user 4.0K Sep 13 07:06 trinity_assemblies
    -rw-r--r-- 1 user user  44M Sep 13 07:08 dm3.fasta.masked_BT2.rev.1.bt2
    -rw-r--r-- 1 user user  30M Sep 13 07:08 dm3.fasta.masked_BT2.rev.2.bt2
    -rw-r--r-- 1 user user  44M Sep 13 07:08 dm3.fasta.masked_BT2.1.bt2
    -rw-r--r-- 1 user user  30M Sep 13 07:08 dm3.fasta.masked_BT2.2.bt2
    -rw-r--r-- 1 user user 505K Sep 13 07:07 dm3.fasta.masked_BT2.3.bt2
    -rw-r--r-- 1 user user  30M Sep 13 07:07 dm3.fasta.masked_BT2.4.bt2
    -rw-r--r-- 1 user user 165M Sep 13 07:07 dm3.fasta.masked
    -rw-r--r-- 1 user user 5.2M Sep 13 07:06 repeatsToMask_coverage-0.5.gff3
    -rw-r--r-- 1 user user 4.6K Sep 13 07:06 allcandidates_coverage-0.5.gff3


The candidate Transposable Elements can be found in the **allcandidates_coverage-0.5.gff3**, those that were removed from the genome in the **repeatsToMask_coverage-0.5.gff3** file. The genome, with the repeats in **repeatsToMask_coverage-0.5.gff3** removed, is in the **dm3.fasta.masked**. Additional files for Bowtie 2 are also generated (**\*\_BT2\*.bt2**).

### CONTACT
_________________________________________________

Please send any inquiries about usage and/or bugs to TEcandidates@gmail.com




### REFERENCES
_________________________________________________

	Ohtani H, Iwasaki YW, Shibuya A, Siomi H, Siomi MC, Saito K. (2013). DmGTSF1 is necessary for Piwi-piRISC-mediated transcriptional transposon silencing in the Drosophila ovary. Genes Dev. 2013 Aug 1;27(15):1656-61. doi: 10.1101/gad.221515.113.



