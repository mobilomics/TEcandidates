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
	-m=*|--mode=*)
        mode="${i#*=}"
        ;;
	-fq=*|--fastq-files=*)
	fastqFiles="${i#*=}"
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
	-N=*)
        N="${i#*=}"
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

if [ -z $mode ]
then
        echo "Please specify a mode: SE for Single-End Reads or PE for Paired-End Reads"
        exit
fi

if [ -z $fastqFiles ]
then
        echo "Please specify the folder containing the FASTQ files"
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
####PIPELINE START


if [ "$mode" == "SE" ];
then
	echo "TEcandidates started in Single-End mode"
	fastqFiles=$fastqFiles"/*fastq"
else
	echo "TEcandidates started in Paired-End mode"
	fastqFiles=$fastqFiles"/*_1.fastq"
fi

readFiles=$(readlink -f $fastqFiles)
echo "The following read files will be used for the pipeline:"
echo $readFiles

Genome=$(readlink -f $Genome)
TEannotation=$(readlink -f $TEannotation)

echo "Genome file: $Genome"
echo "Transposable element annotation: $TEannotation"
echo "Coverage: $minimumCoverage"
echo "Memory: $RAM"
echo "Threads: $threads"
echo "Number of candidates per TE instance: $N"

CPU=$threads

currentworkingdir=$(pwd)

TEcandidatesDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

cd $currentworkingdir

outputdir="candidateTE_analysis_coverage-"$minimumCoverage"_length-"$minimumLength"_N-"$N

if [ ! -d $outputdir ];
then
	echo "Creating $outputdir ..."
	mkdir $outputdir
else
	echo "Output directory already exists"
fi

cd $outputdir


echo ""
echo "#####PRE-MAPPING READS"

BT2_basename=${Genome/.fasta/_BT2}

echo "Building BT2 index for original genome"
bowtie2-build $Genome $BT2_basename -q --threads $threads
echo "Done building BT2 index"

filteredReadFiles=""

for readFile in $readFiles;
do

	readFileBasename=$(basename $readFile)

	if [ "$mode" == "SE" ];
        then
		readFileBasename=$(basename $readFile)

		    readFileBasename=${readFileBasename%.*}
    samFile=${readFileBasename}".sam"
    bowtie2Summary=${readFileBasename}".bt2_summary"
    bamFile=${readFileBasename}".bam"
    filteredBamFile=${readFileBasename}"_filtered.bam"
    filteredFastq=${readFileBasename}"_filtered.fastq"

		if [ -e $filteredFastq ];
		then
			echo "Filtered Fastq for $readFile exists! Filename: $filteredFastq"
		else
			echo "Mapping reads from $readFile to genome"
			bowtie2 -N 1 --threads $threads -x $BT2_basename -U $readFile -S $samFile 2> $bowtie2Summary
			echo "Pre-mapping finished"

			echo "Creating filtered reads file"
			samtools view $samFile --threads $threads -b|bedtools intersect -a - -b $TEannotation -u > $filteredBamFile
			samtools fastq $filteredBamFile --threads $threads > $filteredFastq

			#Conversion to FASTA input (if required)
			#awk '(NR%4==1||NR%4==2){print $0}' $filteredFastq|sed 's/^@/>/g' > temp.txt
			#mv temp.txt $filteredFastq

			rm $samFile $filteredBamFile
			echo "Done creating filtered reads file $filteredFastq"
		fi

                filteredFastq=$(readlink -f $filteredFastq)
                filteredReadFiles=$filteredReadFiles""$filteredFastq" "


        else
		mate1=$readFile
		mate2=${mate1/_1.fastq/_2.fastq}

		mate1Basename=$(basename $mate1)
		samFile=${mate1Basename/_1.fastq/.sam}
		bowtie2Summary=${mate1Basename/_1.fastq/.bt2_summary}

		bamFile=${samFile/.sam/.bam}
		filteredSamFile=${samFile/.sam/_filtered.sam}
		#filteredBamFile=${samFile/.sam/_filtered.bam}
		filteredFastq1=${bamFile/.bam/_filtered_1.fastq}
		filteredFastq2=${bamFile/.bam/_filtered_2.fastq}

		#Make BED from RepeatMasker GFF3
		bedFile=${TEannotation/.gff3/.bed}
		awk 'BEGIN{FS=OFS="\t"}{$4=$4-1;print $1,$4,$5,$9,$6,$7}' $TEannotation > $bedFile


		IDSfile="IDs_that_overlap.txt"
		#threads=100
		#BT2_basename="Xla.v91_BT2"

		if [ -e $filteredFastq1 ];
		then
			echo "Filtered Fastq for $mate1 exists! Filename: $filteredFastq1"  
		else
			#bowtie2 -N 1 --threads $threads -x $BT2_basename -1 $mate1 -2 $mate2 2>$bowtie2Summary|samtools view --threads $threads -b|samtools sort --threads $threads -o $bamFile
			samtools view --threads $threads -L $bedFile $bamFile|awk '{print $1}'|sort -u > $IDSfile
			samtools view -H $bamFile > $filteredSamFile
			samtools view --threads $threads $bamFile|grep -Ff $IDSfile --word-regexp >> $filteredSamFile
			samtools fastq --threads $threads -1 $filteredFastq1 -2 $filteredFastq2 -N $filteredSamFile
			rm $IDSfile $filteredSamFile $bedFile
		fi

		filteredFastq=$(readlink -f $filteredFastq1)
                filteredReadFiles=$filteredReadFiles""$filteredFastq1" "	
	fi        	 

        
done
######
echo $filteredReadFiles

echo "#####GENERATING TRINITY ASSEMBLIES"

if [ ! -d trinity_assemblies ];
then
	mkdir trinity_assemblies
fi

#for readFile in $readFiles;
for readFile in $filteredReadFiles;
do

	echo "Generating Trinity assembly with $readFile"
	basename=$(basename $readFile)
	timeOutput="trinity_assemblies/"${basename/.fastq/.time}
	trinityOutput="trinity_assemblies/"${basename/.fastq/.trinity_assembly.fasta}

	if [ -e $trinityOutput ];
	then
		echo "$trinityOutput exists. Skipping this file"
	else
		if [ "$mode" == "SE" ];
		then
			cmd="/usr/bin/time -f \"%E real\n%U user\n%S sys\n%K memory\" -o $timeOutput Trinity --seqType fq --max_memory $RAM --CPU $CPU --bflyHeapSpaceMax $bflyHeapMax --bflyCPU $CPU --single $readFile --full_cleanup"
			echo "CMD: $cmd"
			echo -e "\n"
			#Trinity --seqType fq --single $readFile --max_memory $RAM --CPU $CPU --bflyHeapSpaceMax $RAM --bflyCPU $CPU --full_cleanup
			/usr/bin/time -f "%E real\n%U user\n%S sys\n%K memory" -o $timeOutput Trinity --seqType fq --max_memory $RAM --CPU $CPU --bflyHeapSpaceMax $bflyHeapMax --bflyCPU $CPU --single $readFile --full_cleanup
		else
			mate1=$readFile
			mate2=${mate1/_1.fastq/_2.fastq}
			cmd="/usr/bin/time -f \"%E real\n%U user\n%S sys\n%K memory\" -o $timeOutput Trinity --seqType fq --max_memory 128G --CPU 12 --bflyHeapSpaceMax 10G --bflyCPU 12 --left $mate1 --right $mate2 --full_cleanup"
			#Trinity --seqType fq --max_memory $RAM --CPU $CPU --bflyHeapSpaceMax $RAM --bflyCPU $CPU --left $mate1 --right $mate2 --full_cleanup
			/usr/bin/time -f "%E real\n%U user\n%S sys\n%K memory" -o $timeOutput Trinity --seqType fq --max_memory $RAM --CPU $CPU --bflyHeapSpaceMax $bflyHeapMax --bflyCPU $CPU --left $mate1 --right $mate2 --full_cleanup
		fi
	
		cmd="mv trinity_out_dir.Trinity.fasta $trinityOutput"
		echo -e "CMD: $cmd\n"
		mv trinity_out_dir.Trinity.fasta $trinityOutput

		if [ $? -ne 0 ];
		then
			"Could not execute $cmd"
			exit
		fi
		echo "Trinity output stored in $trinityOutput"
	fi

	parsed_trinityOutput=$trinityOutput".parsed"
        echo ">CMD: perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' $trinityOutput > $parsed_trinityOutput"
        perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' $trinityOutput > $parsed_trinityOutput


done

###
echo "#####MAPPING ASSEMBLIES INTO GENOME"

#BT2_basename=$(basename $Genome)
#BT2_basename=${BT2_basename/.fasta/_BT2}


cmd="bowtie2-build --threads $threads $Genome $BT2_basename"
echo -e "CMD: $cmd"
#bowtie2-build --threads $threads $Genome $BT2_basename

for assemblyFile in trinity_assemblies/*.trinity_assembly.fasta.parsed;
do
	echo "Mapping $assemblyFile into $BT2_basename"
	samFile=${assemblyFile/fasta.parsed/sam}
	bowtie2Summary=${assemblyFile/fasta.parsed/bt2_summary}
	cmd="bowtie2 -N 1 --threads $threads -f -x $BT2_basename -U $assemblyFile -S $samFile 2> $bowtie2Summary"
	echo -e "CMD: $cmd\n"

	if [ ! -e $bowtie2Summary ]
	then
		bowtie2 -N 1 --threads $threads -f -x $BT2_basename -U $assemblyFile -S $samFile 2> $bowtie2Summary
	else
		if [ $(grep -c "overall" $bowtie2Summary) -eq 0 ]
		then
			bowtie2 -N 1 --threads $threads -f -x $BT2_basename -U $assemblyFile -S $samFile 2> $bowtie2Summary
		else
			echo "$assemblyFile already mapped. Skipping this step"
		fi	
	fi

	bamFile=${samFile%.*}".bam"
	bedFile=${samFile%.*}".bed"
#	bamFile=${samFile/sam/bam}

	if [ -e $samFile ]
	then
		cmd="samtools view -o $bamFile --threads $threads -b $samFile"
		echo -e "CMD: $cmd\n"
		samtools view -o $bamFile --threads $threads -b $samFile
		rm $samFile
	fi

#	bedFile=${samFile/sam/bed}
	cmd="bedtools bamtobed -i $bamFile > $bedFile"
	echo -e "CMD: $cmd\n"
	bedtools bamtobed -i $bamFile > $bedFile

done

cmd="rm $BT2_basename*"
echo -e "CMD: $cmd\n"
rm $BT2_basename*
####

echo "#####GENERATING CANDIDATES"

for trinityBed in trinity_assemblies/*.trinity_assembly.bed;
do
        echo "Generating candidate TE file for $trinityBed..."
        candidateTEfilename=${trinityBed/.trinity_assembly.bed/.candidate_TEs}
        ##candidateTEfilename=$outputdir"/"$candidateTEfilename
        #bedtools coverage -s -a $TEannotation -b $trinityBed |awk -v minimumCoverage=$minimumCoverage '($(NF-1)>900 && $NF>=minimumCoverage) {print $0}' > $candidateTEfilename

	cmd="bedtools coverage -s -a $TEannotation -b $trinityBed |awk -v minimumCoverage=$minimumCoverage -v minimumLength=$minimumLength '(\$(NF-1)>=minimumLength && \$NF>=minimumCoverage) {print \$0}' > $candidateTEfilename"
	echo -e "CMD: $cmd\n"
	bedtools coverage -s -a $TEannotation -b $trinityBed |awk -v minimumCoverage=$minimumCoverage -v minimumLength=$minimumLength '($(NF-1)>=minimumLength && $NF>=minimumCoverage) {print $0}' > $candidateTEfilename
	echo -e "$candidateTEfilename generated!\n\n"
done

cd "trinity_assemblies"
cat *.candidate_TEs > allcandidates.gff3
sort -k9,9d -k12,12n -k13,13n allcandidates.gff3 > allcandidates.gff3.sorted
finalfile="allcandidates_coverage-"$minimumCoverage"_length-"$minimumLength"_N-"$N".gff3"

echo $finalfile

#################

allcandidatesFile="allcandidates.gff3.sorted"
idsFile="uniqueIDS.txt"
groupedcandidatesFile=$allcandidatesFile".grouped"

bedtools groupby -g 1,2,3,4,5,6,7,8,9 -c 10,11,12,13 -o last -i $allcandidatesFile > $groupedcandidatesFile

awk '{print $9}' $allcandidatesFile|sort -u > $idsFile
while read line
do
        pattern=$line
        cmd="grep -w $pattern $groupedcandidatesFile|tail -n $N >> $finalfile"
        #echo "$cmd"
        grep -w $pattern $groupedcandidatesFile|tail -n $N >> $finalfile
done < "$idsFile"

#################
#bedtools groupby -g 9 -c 1,2,3,4,5,6,7,8,9 -o last -i allcandidates.gff3.sorted|awk 'BEGIN{OFS="\t"} {$10=$10"/ChrID="$2$8; print $2,$3,$4,$5,$6,$7,$8,$9,$10}' > $finalfile

awk 'BEGIN{OFS="\t"} {$9=$9"/ChrID="$1$7"_"$4"-"$5; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' $finalfile > temp.txt
#awk 'BEGIN{OFS="\t"} {$9=$9"_"$4"-"$5; print $1,$2,$3,$4,$5,$6,$7,$8,$9}' $finalfile > temp.txt
mv temp.txt $finalfile 

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

