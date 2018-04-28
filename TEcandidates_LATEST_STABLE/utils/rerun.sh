#!/bin/bash

for i in "$@"
do
	case $i in
	-g=*|--genome=*)
	Genome="${i#*=}"
	;;
	-t=*|--threads=*)
	threads="${i#*=}"
	;;
        -c=*|--coverage=*)
        minimumCoverage="${i#*=}"
        ;;
	-l=*|--length=*)
        minimumLength="${i#*=}"
        ;;
	-te=*|--te-annotation=*)
        TEannotation="${i#*=}"
        ;;
	-r=*|--RAM=*)
        RAM="${i#*=}"
        ;;
	-d=*|--directory=*)
        previousRunDirectory="${i#*=}"
        ;;
	esac
done

####DEPENDENCIES VALIDATION


echo "Checking if Bedtools v2.25 is installed.."
if hash bedtools 2>/dev/null; then
        bedtools_version=$(bedtools --version)
        if [[ $bedtools_version == "bedtools v2.25.0" ]]; then
                echo $bedtools_version
                echo "Bedtools v2.25 is installed!"
        fi
else
        echo "Bedtools v2.25 was not found. Please install it and try again."
        exit
fi
echo ""


echo "Checking if Bioperl is installed.."
bioperl_version=$(perl -MBio::Root::Version -e 'print $Bio::Root::Version::VERSION,"\n"')
if [[ $bioperl_version != *"Can't locate"* ]]; then
        echo $bioperl_version
        echo "Bioperl is installed!"
else
        echo "Bioperl was not found. Please install it and try again."
        exit
fi
echo ""


echo "Checking if Bowtie v2.3 is installed.."
if hash bowtie2 2>/dev/null; then
        bowtie2_version=$(bowtie2 --version|head -n 1)
        if [[ $bowtie2_version == *"bowtie2-align-s version 2.3.0"* ]]; then
                echo $bowtie2_version
                echo "Bowtie v2.3 is installed!"
        fi
else
        echo "Bowtie2 v2.3 was not found. Please install it and try again."
        exit
fi
echo ""


####INPUT VALIDATION

if [ -z $Genome ]
then
	echo "Please specify a Genome FASTA file"
	exit
fi

if [ -z $threads ]
then
        echo "Please specify a number of threads to use during the pipeline execution"
        exit
fi

if [ -z $RAM ]
then
        echo "Please specify the amount of RAM (in GB) to use during the pipeline execution"
        exit
fi

if [ -z $previousRunDirectory ]
then
        echo "Please specify the folder containing the previous run results"
        exit
fi

if [ -z $minimumCoverage ]
then
        echo "Please specify the coverage at which candidate TEs will be selected"
        exit
fi

if [ -z $minimumLength ]
then
        echo "Please specify the length at which candidate TEs will be selected"
        exit
fi

if [ -z $TEannotation ]
then
        echo "Please specify a TE annotation file in GFF3 format"
        exit
fi


CPU=$threads
bflyHeapMax=$(($RAM/$CPU))
bflyHeapMax=$bflyHeapMax"G"
RAM=$RAM"G"
previousRunDirectory=$(readlink -f $previousRunDirectory)

####PIPELINE START

Genome=$(readlink -f $Genome)
TEannotation=$(readlink -f $TEannotation)

echo "Genome file: $Genome"
echo "Transposable element annotation: $TEannotation"
echo "Coverage: $minimumCoverage"
echo "Memory: $RAM"
echo "Threads: $threads"
echo "Previous run directory: $previousRunDirectory"

CPU=$threads

currentworkingdir=$(pwd)

TEcandidatesDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $currentworkingdir

outputdir="candidateTE_analysis_coverage-"$minimumCoverage"_length-"$minimumLength

if [ ! -d $outputdir ];
then
	echo "Creating $outputdir ..."
	mkdir $outputdir
else
	echo "Output directory already exists"
fi

outputdir=$(readlink -f $outputdir)
cd $outputdir
mkdir trinity_assemblies

echo "#####GENERATING CANDIDATES"

for trinityBed in $previousRunDirectory/trinity_assemblies/*.trinity_assembly.bed;
do
        echo "Generating candidate TE file for $trinityBed..."
	
        candidateTEfilename=${trinityBed/.trinity_assembly.bed/.candidate_TEs}
	candidateTEfilename=$(basename $candidateTEfilename)
	candidateTEfilename=$outputdir"/trinity_assemblies/"$candidateTEfilename	

	
	cmd="bedtools coverage -s -a $TEannotation -b $trinityBed |awk -v minimumCoverage=$minimumCoverage -v minimumLength=$minimumLength '(\$(NF-1)>=minimumLength && \$NF>=minimumCoverage) {print \$0}' > $candidateTEfilename"
	echo -e "CMD: $cmd\n"
	bedtools coverage -s -a $TEannotation -b $trinityBed |awk -v minimumCoverage=$minimumCoverage -v minimumLength=$minimumLength '($(NF-1)>=minimumLength && $NF>=minimumCoverage) {print $0}' > $candidateTEfilename

	echo -e "$candidateTEfilename generated!\n\n"
done

cd "trinity_assemblies"
cat *.candidate_TEs > allcandidates.gff3

#Sort the file with the coverage info, according to their ID (column 9), then by the TE instance length (column 12) and by the coverage by a de novo contig (column 13). 
sort -k9,9d -k12,12n -k13,13n allcandidates.gff3 > allcandidates.gff3.sorted
finalfile="allcandidates_coverage-"$minimumCoverage"_length-"$minimumLength".gff3"

echo $finalfile

#Column 9 is where the ID of a TE is. Here, the files with the coverages per TE instance is grouped by their ID, in order to find the instance of a given TE that has the greates coverage. Then, the location is appended to the ID (this might be redundant, as this is a GFF3 file after all).
bedtools groupby -g 9 -c 1,2,3,4,5,6,7,8,9 -o last -i allcandidates.gff3.sorted|awk 'BEGIN{OFS="\t"} {$10=$10"/ChrID="$2$8; print $2,$3,$4,$5,$6,$7,$8,$9,$10}' > $finalfile

mv $finalfile ../
cd ..

echo "#####GENERATING MASKED GENOME AND FINAL FILES"
repeatsToMask="repeatsToMask_coverage-"$minimumCoverage"_length-"$minimumLength".gff3"
cmd="bedtools intersect -v -a $TEannotation -b $finalfile > $repeatsToMask"
echo -e "CMD: $cmd\n"
bedtools intersect -v -a $TEannotation -b $finalfile > $repeatsToMask

maskedGenome=$(basename $Genome)
maskedGenome=$maskedGenome".masked"
BT2_basename=$maskedGenome"_BT2"
cmd="bedtools maskfasta -fi $Genome -fo $maskedGenome -bed $repeatsToMask -mc X"
echo -e "CMD: $cmd\n"
bedtools maskfasta -fi $Genome -fo $maskedGenome -bed $repeatsToMask -mc X

cmd="bowtie2-build --threads $threads $maskedGenome $BT2_basename"
echo -e "CMD: $cmd\n"
bowtie2-build --threads $threads $maskedGenome $BT2_basename

