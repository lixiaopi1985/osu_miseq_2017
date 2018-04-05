#!/bin/bash

# this  shell script integrates following steps:
# 1 step: fastqc check all the reads quality
# (assuming they are all in one folder).
# 2 step: trim the reads with trimmomatic
# 3 step: check reads with fastqc again (generates another dir)
# 4 step: group the reads into different folders according to
# the sample they belong to
# 5 step: assemble the pair ended reads for each sample
# 6 step: combine reads

usage() { echo "Usage: $0 [-h, --help] [ -p, --path: path file name] [ -o, --output: the folder name to trimmed reads] [-s, --sample: name for the file that contains sample names]" >&2; exit 1;}

while getopts ":p:o:s:r:" args; do
	case "${args}" in
		p|--path)
			paths=${OPTARG}
			;;
		o|--output)
			outfolder=${OPTARG}
			;;
		s|--sample)
			sample=${OPTARG}
			;;
		r|--primer)
			primer=${OPTARG}
			;;
		h|--help)
			usage
			;;
		*)
			usage
			;;
	esac
done
shift $((OPTIND-1))


printf "Running step 1: Checking reads quality ... \n"


# Find raw_reads
# Find illumina seq adapters
path=$( grep "raw_reads" $paths | cut -d ':' -f 2)
adapters_path=$( grep "adapters" $paths | cut -d ':' -f 2)


printf "\n"
mkdir step1_rawRds_fastqc_report

cd step1_rawRds_fastqc_report


printf "\n"


echo "Combining all reads for fastqc, it will take a few minutes"


cat $path/*.gz > all_raw_reads.fastq.gz & P1=$!

spin='-\|/'

i=0
while kill -0 $P1 2> /dev/null
do
	i=$(( (i+1) %4 ))
	printf "\r${spin:$i:1}"
	sleep .1
done

### setting wait to prevent other commands run before this one finishes
wait $P1

printf "\n"

fastqc all_raw_reads.fastq.gz & P2=$!

spin='-\|/'

i=0
while kill -0 $P2 2> /dev/null
do
	i=$(( (i+1) %4 ))
	printf "\r${spin:$i:1}"
	sleep .1
done

wait $P2

printf "\nRemoving all_raw_reads.fastq.gz to save disk space\n"
rm all_raw_reads.fastq.gz
cd ..

printf "\nFastqc on the raw reads completed!\n\n"

#############################################################################
################ trimming ###################################################
#############################################################################
#############################################################################

printf "Step 2: use trimmomatic to trim the raw reads ...\n"


# entering trim reads folder
mkdir $outfolder
printf "\n"
printf "Setting parameters for trimmomatic\n"

echo -n "LEADING >>"
read leading
echo -n "TRAILING >>"
read trailing
echo -n "SLIDINGWINDOW WIDTH >>"
read S_width
echo -n "SLIDINGWINDOW QUALITY >>"
read S_q
echo -n "MINLEN >>"
read minlen
echo -n "Thread >>"
read threads

printf "\nStarting trimming\n"

# number of fq in raw reads file
nReads=$( ls $path | wc -l )

for (( i=1;i<=$nReads/2;i++))
do
	label=S$i
	printf "Working on Sample $label \n"
	# right read
	right=$( ls $path | grep .*_${label}_.*_R1_.* )
	# left read
	left=$( ls $path | grep .*_${label}_.*_R2_.* )

	java -jar $TRIMO PE $path/$right $path/$left \
$outfolder/"trimmed_paired_${right%%.gz}" \
$outfolder/"trimmed_unpaired_${right%%.gz}" \
$outfolder/"trimmed_paired_${left%%.gz}" \
$outfolder/"trimmed_upaired_${left%%.gz}" \
-threads $threads \
ILLUMINACLIP:$adapters_path:2:30:10:8:true \
LEADING:$leading \
TRAILING:$trailing \
SLIDINGWINDOW:$S_width:$S_q \
MINLEN:$minlen

done

printf "\nTrimming Completed!\n"


#########################################################################
################## Check trimmed quality ################################
#########################################################################


printf "Now Fastqc is checking trimmed reads quality ... \n"


mkdir step3_check_trimmed_QC

cd step3_check_trimmed_QC

printf "\n"


echo "Combining all paired reads for fastqc, it will take a few minutes"
cat ../$outfolder/*paired*.fastq > all_paired_trimmed_reads.fastq & P1=$!

spin='-\|/'

i=0
while kill -0 $P1 2> /dev/null
do
	i=$(( (i+1) %4 ))
	printf "\r${spin:$i:1}"
	sleep .1
done

### setting wait to prevent other commands run before this one finishes
wait $P1

printf "\n"

echo "Combining all unpaired reads for fastqc, it will take a few minutes"
cat ../$outfolder/*unpaired*.fastq > all_unpaired_trimmed_reads.fastq & P2=$!
spin='-\||'
i=0
while kill -0 $P2 2> /dev/null
do
	i=$(( (i+1) %4 ))
	printf "\r${spin:$i:1}"
	sleep .1
done
wait $P2

printf "Fastqc on all paired read\n"
fastqc all_paired_trimmed_reads.fastq & P3=$!
wait $P3
# multiqc report

multiqc .

printf "Fastqc on all unpaired reads\n"
fastqc all_unpaired_trimmed_reads.fastq

multiqc .

printf "\nRemoving fastq to save disk space\n"
rm all_paired_trimmed_reads.fastq
rm all_unpaired_trimmed_reads.fastq

cd ..

############################################################################
############################## Group reads #################################
############################################################################


printf "Grouping reads according to barcodes ... \n"

echo "\n"
echo "\n"

if [ -e "step4_trimread_grouping" ];then
	:
else
	mkdir step4_trimread_grouping
fi

cd step4_trimread_grouping

step4_p=$( pwd )


## sample file needs to contain D1_32, V1_32, K1_16, K17_32
while IFS='' read -r line || [[ -n "$line" ]];do
	folder=$line
	if [ -e $folder ];then
		:
	else
		mkdir $folder
		mkdir unpaired_$folder
	fi
done < ../$sample

cd ..


#get group directories path
v1_32_p=$step4_p/$( ls $step4_p | grep "^V" )
d1_32_p=$step4_p/$( ls $step4_p | grep "^D"  )
k1_16_p=$step4_p/$( ls $step4_p | grep "^K1_" )
k17_32_p=$step4_p/$( ls $step4_p | grep "^K17_" )



# get group unpaired
v1_32_p_un=$step4_p/$( ls $step4_p | grep "unpaired_V" )
d1_32_p_un=$step4_p/$( ls $step4_p | grep "unpaired_D"  )
k1_16_p_un=$step4_p/$( ls $step4_p | grep "unpaired_K1_" )
k17_32_p_un=$step4_p/$( ls $step4_p | grep "unpaired_K17_" )



cd $step4_p


for i in $(ls ../$outfolder)
do
        if [[ "$i" =~ .*_paired_.*D-[0-9]{1,2} ]];then
                mv ../$outfolder/$i $v1_32_p

        elif [[ "$i" =~ .*_unpaired_.*D-[0-9]{1,2} ]];then
                mv ../$outfolder/$i $v1_32_p_un

        elif [[ "$i" =~ .*_paired_.*K-([1-9]|1[0-6])_ ]];then
                mv ../$outfolder/$i $k1_16_p

        elif [[ "$i" =~ .*_unpaired_.*K-([1-9]|1[0-6])_ ]];then
                mv ../$outfolder/$i $k1_16_p_un

        elif [[ "$i" =~ .*_paired_.*K-(1[7-9]|[2-9][0-9])_ ]];then
                mv ../$outfolder/$i $k17_32

        elif [[ "$i" =~ .*_unpaired_.*K-(0[7-9]|[2-9][0-9])_ ]];then
                mv ../$outfolder/$i $k17_32_un

        elif [[ "$i" =~ .*_paired_.*V-[0-9]{1,2} ]];then
                mv ../$outfolder/$i $v1_32_p

        elif [[ "$i" =~ .*_unpaired_.*V-[0-9]{1,2} ]];then
                mv ../$outfolder/$i $v1_32_p_un

        fi
done


cd ..





# primers
trnL_f=$( grep "^trnL_f" $primer | cut -d ';' -f 2 )
UAA_h_comp=$( grep "^UAA_h_comp" $primer | cut -d ';' -f 2 )
UAA_h=$( grep "^UAA_h;" $primer | cut -d ';' -f 2 )
ITS86_f=$( grep "^ITS86_f" $primer | cut -d ';' -f 2 )
ITS4_r_comp=$( grep "^ITS4_r_comp" $primer | cut -d ';' -f 2 )
ITS4_r=$( grep "^ITS4_r;" $primer | cut -d ';' -f 2 )
ZBJ_f=$( grep "^ZBJ_f" $primer | cut -d ';' -f 2 )
ZBJ_r_comp=$( grep "^ZBJ_r_comp" $primer | cut -d ';' -f 2 )
ZBJ_r=$( grep "^ZBJ_r;" $primer | cut -d ';' -f 2 )




#############################################################################
######################### assembly ##########################################
#############################################################################

if [ ! -e "step5_assembly" ];then
	mkdir step5_assembly
fi

cd step5_assembly
assembly_p=$( pwd )


## make sub folders
while IFS='' read -r line || [[ -n "$line" ]];do
        folder=$line
	if [ ! -e "$folder" ];then
        	mkdir $folder
		mkdir fasta_${folder}
	fi
done < ../$sample

# exit step5
cd ..

printf "Assemble ... \n"

v1_32_p2=$assembly_p/$( ls $assembly_p | grep "^V" )
d1_32_p2=$assembly_p/$( ls $assembly_p | grep "^D"  )
k1_16_p2=$assembly_p/$( ls $assembly_p | grep "^K1_" )
k17_32_p2=$assembly_p/$( ls $assembly_p | grep "^K17_" )


v1_32_p2_fa=$assembly_p/$( ls $assembly_p | grep "^fasta_V" )
d1_32_p2_fa=$assembly_p/$( ls $assembly_p | grep "^fasta_D"  )
k1_16_p2_fa=$assembly_p/$( ls $assembly_p | grep "^fasta_K1_" )
k17_32_p2_fa=$assembly_p/$( ls $assembly_p | grep "^fasta_K17_" )



# go to step4 look for v samples


cd $v1_32_p

for i in {1..32}
do
        Fw=$(ls | grep ^trimmed_.*-V-${i}_.*_R1_.* )
        Rv=$(ls | grep ^trimmed_.*-V-${i}_.*_R2_.* )

        printf "Applying cope basic to assemble the pair ends, sample V${i}...\n"
        cope -a $Fw -b $Rv -o $v1_32_p2/V_${i}_coped.fastq -2 $v1_32_p2/unAssembled_$Fw -3 $v1_32_p2/unAssembled_$Rv -m 0 -u 310 -l 10 -c 0.75 >$v1_32_p2/cope_V$i.log 2>$v1_32_p2/cope_V$i.error


        printf "Trimming off trnL and UAA primers in the sequence...\n"
        cutadapt -a $trnL_f...$UAA_h_comp $v1_32_p2/V_${i}_coped.fastq > $v1_32_p2/cutadapt1_V_${i}_coped.fastq
        cutadapt -g $trnL_f -a $UAA_h_comp $v1_32_p2/cutadapt1_V_${i}_coped.fastq > $v1_32_p2/cutadapt2_V_${i}_coped.fastq
        cutadapt -g $UAA_h $v1_32_p2/cutadapt2_V_${i}_coped.fastq > $v1_32_p2/cutadapt3_V_${i}_coped.fastq


done


printf "assemble D1_32 paired end reads...\n"
cd $d1_32_p

for i in {1..32}
do
        Fw=$(ls | grep ^trimmed_.*-D-${i}_.*_R1_.* )
        Rv=$(ls | grep ^trimmed_.*-D-${i}_.*_R2_.* )


        printf "Applying cope basic to assemble the pair ends, sample D${i}...\n"
        cope -a $Fw -b $Rv -o $d1_32_p2/D_${i}_coped.fastq -2 $d1_32_p2/unAssembled_$Fw -3 $d1_32_p2/unAssembled_$Rv -m 0 \
-u 310 -l 10 -c 0.75 >$d1_32_p2/cope_D$i.log 2>$d1_32_p2/cope_D$i.error
        printf "Trimming off zbj_f and zbj_r primers in the sequence...\n"
        cutadapt -a $ZBJ_f...$ZBJ_r_comp $d1_32_p2/D_${i}_coped.fastq > $d1_32_p2/cutadapt1_D_${i}_coped.fastq
        cutadapt -g $ZBJ_f -a $ZBJ_r_comp $d1_32_p2/cutadapt1_D_${i}_coped.fastq > $d1_32_p2/cutadapt2_D_${i}_coped.fastq
        cutadapt -g $ZBJ_r $d1_32_p2/cutadapt2_D_${i}_coped.fastq > $d1_32_p2/cutadapt3_D_${i}_coped.fastq

done



printf "assemble K1_16 paired end reads...\n"
cd $k1_16_p

for i in {1..16}
do
        Fw=$(ls | grep ^trimmed_.*-K-${i}_.*_R1_.* )
        Rv=$(ls | grep ^trimmed_.*-K-${i}_.*_R2_.* )


        printf "Applying cope basic to assemble the pair ends, sample K${i}...\n"
        cope -a $Fw -b $Rv -o $k1_16_p2/K_${i}_coped.fastq -2 $k1_16_p2/unAssembled_$Fw -3 $k1_16_p2/unAssembled_$Rv -m 0 \
-u 310 -l 10 -c 0.75 >$k1_16_p2/cope_K$i.log 2>$k1_16_p2/cope_K$i.error
        printf "Trimming off trnL and UAA primers in the sequence...\n"
        cutadapt -a $trnL_f...$UAA_h_comp $k1_16_p2/K_${i}_coped.fastq > $k1_16_p2/cutadapt1_K_${i}_coped.fastq
        cutadapt -g $trnL_f -a $UAA_h_comp $k1_16_p2/cutadapt1_K_${i}_coped.fastq > $k1_16_p2/cutadapt2_K_${i}_coped.fastq
        cutadapt -g $UAA_h $k1_16_p2/cutadapt2_K_${i}_coped.fastq > $k1_16_p2/cutadapt3_K_${i}_coped.fastq

done


printf "assemble K17_32 paired end reads...\n"

cd $k17_32_p

for i in {17..32}
do
        Fw=$(ls | grep ^trimmed_.*-K-${i}_.*_R1_.* )
        Rv=$(ls | grep ^trimmed_.*-K-${i}_.*_R2_.* )


        printf "Applying cope basic to assemble the pair ends, sample K${i}...\n"
        cope -a $Fw -b $Rv -o $k17_32_p2/K_${i}_coped.fastq -2 $k17_32_p2/unAssembled_$Fw -3 $k17_32_p2/unAssembled_$Rv -m 0 \
-u 310 -l 10 -c 0.75 >$k17_32_p2/cope_K$i.log 2>$k17_32_p2/cope_K$i.error
        printf "Trimming off ITS86F and ITS4 primers in the sequence...\n"
        cutadapt -a $ITS86_f...$ITS4_r_comp $k17_32_p2/K_${i}_coped.fastq > $k17_32_p2/cutadapt1_K_${i}_coped.fastq
        cutadapt -g $ITS86_f -a $ITS4_r_comp $k17_32_p2/cutadapt1_K_${i}_coped.fastq > $k17_32_p2/cutadapt2_K_${i}_coped.fastq
        cutadapt -g $ITS4_r $k17_32_p2/cutadapt2_K_${i}_coped.fastq > $k17_32_p2/cutadapt3_K_${i}_coped.fastq
done

cd ../../


############################################################################################################################################  step 6 ######################################################################################################################################################################


printf "Converting Fastq file to fasta file\n"

cd $v1_32_p2

echo "In the directory >>>> $v1_32_p2\n" >&2


for i in {1..32};do
	FQ=$( ls | grep ^cutadapt3_V_${i}_coped.fastq )
	echo $FQ
	python3 ../fastq2fasta.py -i $FQ -o $v1_32_p2_fa -d y
done



cd $d1_32_p2
echo "In the directory >>>> $d1_32_p2\n" >&2
for i in {1..32};do
        FQ=$( ls | grep ^cutadapt3_D_${i}_coped.fastq )
        echo $FQ
        python3 ../fastq2fasta.py -i $FQ -o $d1_32_p2_fa -d y
done

cd $k1_16_p2
echo "In the directory >>>> $k1_16_p2\n" >&2
for i in {1..16};do
        FQ=$( ls | grep ^cutadapt3_K_${i}_coped.fastq )
        echo $FQ
        python3 ../fastq2fasta.py -i $FQ -o $k1_16_p2_fa -d y
done



cd $k17_32_p2
echo "In the directory >>>> $k17_32_p2\n" >&2
for i in {17..32};do
        FQ=$( ls | grep ^cutadapt3_K_${i}_coped.fastq )
        echo $FQ
        python3 ../fastq2fasta.py -i $FQ -o $k17_32_p2_fa -d y
done


cd ..

mkdir combined_reads
cd combined_reads

printf "Concatenating all V samples\n"

cat $v1_32_p2_fa/cutadapt3_*.fasta > combined_cutadapt_V1_32.fasta

printf "Concatenating all D samples\n"  

cat $d1_32_p2_fa/cutadapt3_*.fasta > combined_cutadapt_D1_32.fasta

printf "Concatenating all K samples\n"

cat $k1_16_p2_fa/cutadapt3_*.fasta > combined_cutadapt_K1_16.fasta

cat $k17_32_p2_fa/cutadapt3_*.fasta > combined_cutadapt_K17_32.fasta



printf "Combining task has completed!\n"
printf "You can move on to the OTU picking module\n"


 





















