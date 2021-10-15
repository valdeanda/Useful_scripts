#!/bin/bash
# Ian Rambo - July 22, 2021
# Last updated - Oct 14, 2021
# Thirteen... that's a mighty unlucky number... for somebody!


#Purpose: check BAM files for integrity and coordinate sorting

#Usage: bash bam_check.sh -b <directory with BAM files> -o <output directory> -j <number of BAMs to process in parallel> -n "filename pattern for find command"

## Requirements:
## samtools >= 1.9
## GNU parallel

has_command () {
    #Make sure you can execute certain commands
    command -v "$1" >/dev/null 2>&1 || { echo "Requires $1. Ensure that $1 is in your \$PATH or environment."; exit 1; }
}

has_command samtools && has_command parallel
#=============================================================================
## Command-line options
usage="
$(basename "$0"): check BAM files for integrity and coordinate sorting.

where:

   -h --- show this help message
   -b --- parent directory containing BAM files (ending in .bam)
   -o --- output directory for reports and joblogs
   -j --- number of parallel jobs
   -n --- filename pattern to find BAM files to process, e.g. Meg22*.bam

    "

while getopts ':hb:o:j:n:' option; do
    case "${option}" in
    h) echo "$usage"
       exit ;;
    b) bam_dir=${OPTARG};;
    o) out_dir=${OPTARG};;
    j) njobs=${OPTARG};;
    n) bam_name=${OPTARG};;

    :) printf "missing argument for -%s\n" "$OPTARG" >&2
       echo "$usage" >&2
       exit 1 ;;
    \?) printf "illegal option: -%s\n" "$OPTARG" >&2
        echo "$usage" >&2
        exit 1 ;;
    esac
done

shift $((OPTIND - 1))
#=============================================================================
# Output joblogs and reports
joblog_dir=${out_dir}/joblogs
report_dir=${out_dir}/reports

test -d $joblog_dir || mkdir -p $joblog_dir
test -d $report_dir || mkdir $report_dir

# Output files showing BAMs with issues
joblog_quickcheck=${joblog_dir}/quickcheck.joblog
report_quickcheck=${report_dir}/report_quickcheck.txt

joblog_sorted=${joblog_dir}/sorted.joblog
report_sorted=${report_dir}/report_sorted.txt

# Logfile
progress_log=${joblog_dir}/progress.log

# Count the number of bam files that will be processed
nbam=0
# BAM files matching a specified pattern
if [ ! -z "$bam_name" ]; then
    nbam=$(find $bam_dir -type f -name "$bam_name" | wc -l)

else
    bam_name="*.bam"
    nbam=$(find $bam_dir -type f -name "$bam_name" | wc -l)
    echo "no filename pattern specified with -n option, using $bam_name -- $(date +"%D %T")" >> $progress_log
fi

#If no BAMs are found, exit
if [ "$nbam" -eq "0" ]; then
    echo "ERROR: $nbam BAM files found in $bam_dir matching pattern $bam_name, exiting -- $(date +"%D %T")" >> $progress_log
    echo "ERROR: $nbam BAM files found in $bam_dir matching pattern $bam_name, exiting..." && exit 1
else
    echo "$njobs parallel tests will be run for $nbam BAM files in $bam_dir matching pattern $bam_name"
    echo "$njobs parallel tests will be run for $nbam BAM files in $bam_dir matching pattern $bam_name -- $(date +"%D %T")" >> $progress_log
fi
#=============================================================================
## Samtools quickcheck - check for BAM integrity
echo "running samtools quickcheck to verify BAM integrity -- $(date +"%D %T")" >> $progress_log

# find $bam_dir -type f -name "*.bam" | \
find $bam_dir -type f -name $bam_name | \
    parallel --joblog $joblog_quickcheck --jobs $njobs samtools quickcheck {} && \
    tail -n +2 $joblog_quickcheck | awk '$7 != 0' > $report_quickcheck && \
    echo "BAM integrity samtools quickcheck test completed -- $(date +"%D %T")" >> $progress_log
#=============================================================================
## Samtools stats - check for proper coordinate sorting
echo "running $njobs parallel samtools stats jobs to check if BAMs are properly sorted -- $(date +"%D %T")" >> $progress_log
#find . -type f -name "*.bam" | \
find $bam_dir -type f -name $bam_name | \
    parallel --joblog $joblog_sorted --jobs $njobs samtools stats {} '|' grep -E \"^SN\" '|' cut -f 2- '|' grep -q -E \"is sorted\:[[:space:]]+1\" && \
    tail -n +2 $joblog_sorted | awk '$7 != 0' > $report_sorted && \
    echo "BAM sorting samtools stats test completed -- $(date +"%D %T")" >> $progress_log

echo "COMPLETED -- $(date +"%D %T")" >> $progress_log
