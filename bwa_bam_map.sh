#!/bin/bash
# Ian Rambo - July 22, 2021
# Thirteen... that's a mighty unlucky number... for somebody!

## Requirements:
## samtools >= 1.9
## bwa >= 0.7.17

## FASTQ must be interleaved

has_command () {
    #Make sure you can execute certain commands
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH or environment."; exit 1; }
}


usage="
$(basename "$0"): map reads to an assembly with bwa mem and samtools. Returns an indexed, coordinate-sorted BAM.

where:

   -h --- show this help message
   -i --- input assembly
   -r --- fastq reads file
   -o --- output directory
   -s --- identifier for output BAM, e.g. assembly_M22-reads_M40
   -b --- threads for bwa mem
   -k --- threads for samtools sort
   -t --- temporary file directory
   -e --- error file directory

    "

while getopts ':hi:r:o:s:b:k:t:e:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    i) assembly=${OPTARG};;
    r) fastq=${OPTARG};;
    o) outdir=${OPTARG};;
    s) sampid=${OPTARG};;
    b) bwa_thread=${OPTARG};;
    k) sam_thread=${OPTARG};;
    t) tmp_dir=${OPTARG};;
    e) error_dir=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))

assembly_basename=$(basename $assembly)
fastq_basename=$(basename $fastq)

has_command bwa && has_command samtools

test -d $outdir || mkdir -p $outdir
test -d $tmp_dir || mkdir -p $tmp_dir
test -d $error_dir || mkdir -p $error_dir

#check if the output identifier was specified
if [ -z "$sampid" ]; then
    echo "no output identifier specified with -s, using assembly and fastq reads basenames"
    sampid=${assembly_basename}_${fastq_basename}
else
    echo "using output identifier $sampid"
fi

#Path for output coordinate sorted BAM
bamfile=${outdir}/${sampid}.coordsort.bam

logfile=${outdir}/logfile.txt

echo "mapping $sampid to $assembly -- $(date +"%D %T")" >> $logfile


time bwa mem -t $bwa_thread -p $assembly $fastq | \
    samtools view -hu -F4 - | \
    samtools sort -O bam -T ${tmp_dir}/${sampid} -@ $sam_thread -o - | \
    tee $bamfile | \
    samtools index - ${bamfile}.bai 2>${error_dir}/${sampid}.stderror && \
    echo "mapping $sampid to $assembly COMPLETED -- $(date +"%D %T")" >> $logfile
