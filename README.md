# Useful scripts

This repository contains useful tips/tricks/scripts that I have picked up over the years. 
The are scripts written in 
- bash 
- awk
- perl 
- R
- python3

Most of the python  scripts were written by my collaborator [Cesar Poot](https://github.com/acpooth) and modify it by me. 


## Dependencies

Before using any of the scripts make sure that you have install the following dependencies

```bash
sudo apt-get install python3 python3-pip python3-matplotlib \
ipython3-notebook python3-mpltoolkits.basemap
sudo pip3 install -U pip
sudo -H pip3 install --upgrade pandas numpy scipy seaborn 
sudo -H pip3 install -U scikit-learn

```

---
# Visualization

## Bubble plot

Scrip to create a bubble chart from any dataframe contanining either normalized or absolute values. 

<img src="https://valdeanda.github.io/Useful_scripts/data_bubbleplot.tab_bubbleplot.png" width="400" height="300" align="right">

```bash
usage: bubble_chart.py [-h] [-im_format {png,pdf,ps,eps,svg,tif,jpg}]
                       [--im_res dpi]
                       filename

positional arguments:
  filename              Input file dataframe i.e abundances profile

optional arguments:
  -h, --help            show this help message and exit
  -im_format {png,pdf,ps,eps,svg,tif,jpg}, -f {png,pdf,ps,eps,svg,tif,jpg}
                        Output format for images [png].
  --im_res dpi, -r dpi  Output resolution for images in dot per inch (dpi)
                        [dpi].
 ```
 
 **Running Bubble plot with example data**
 
 ```bash
  python3 bubble_chart.py  data_bubbleplot.tab -f  png -r 300
 ```
 
**Customize your script**

```python
sns.set(font_scale=1) #change font size  
sns.set_style("whitegrid") #whitegrid to change background to white
plt.figure(figsize=(21,12)) #inches, modify to widen (x) or lengthen (y) --> original was 21,12
plt.tight_layout() #keeps axes names in same figure
bubble_super_mega_and_simpe_plot(df, 20, cmap='bone_r', ylabel='Tax Group (# of genomes)', xlabel='Genes',alpha=0.05)
# alpha = transparency
#cmap= color palete see below
 
#Recomended colors 

cmap='bone_r'
cmap='plasma'
cmap='coolwarm_r'
 
#cmap python = see https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html
```
---

# Heatmap 

<img src="https://valdeanda.github.io/Useful_scripts/clustermap.png" width="300" height="450" align="right">


```bash
python3 heatmap.py -h

usage: heatmap.py [-h] [-f {png,pdf,ps,eps,svg,tif,jpg}] [-r dpi] filename

 This script create a cluster map

positional arguments:
  filename              Input file derived mebs_output with classification

optional arguments:
  -h, --help            show this help message and exit
  -f {png,pdf,ps,eps,svg,tif,jpg}, --im_format {png,pdf,ps,eps,svg,tif,jpg}
                        Output format for images [png].
  -r dpi, --im_res dpi  Output resolution for images in dot per inch (dpi)
                        [dpi].

Example:
    $  python3 heatmap.py data.heatmap.tsv
```

---
## Horizontal Barplot 

<img src="https://valdeanda.github.io/Useful_scripts/barplot.png" width="400" height="350" align="right">

R script that was used to plot the number of achaeal genomes by taxonomy described in [Baker et al., 2020](https://www.nature.com/articles/s41564-020-0715-z)

The input data data_barplot.tab, looks like this, the scripts keeps the specific order that you want your data to sorted.
Dr. Craig Connolly provided valuable input to generate the script. 

```
Phylum	Superphylum	Number of genomes
Heimdallarchaeota	Asgard	66
Lokiarchaeota 	Asgard	63
Unclassified_Asgard	Asgard	37
Thorarchaeota	Asgard	32
Helarchaeota	Asgard	19
```

```R
barplot.R
```


---

## Replace names from a phylogenetic tree 

```bash
perl Replace_tree_names.pl mapping_file tree > renamed_tree
```
---

#  Fasta file processing 

## Split fasta 

Requires biopython

```bash
pip3 install biopython
```

Script that is useful if you have a large fasta file and you want to split it into small files of the same size 

```python
python3 split_fasta.py

usage: split_fasta.py [-h] [-p PARTS] fastafile

Split a fasta file according in almost equal parts based on total base/residue
count. Stores a numpy array that contains the lengths of the sequences in the
file

positional arguments:
  fastafile             Fasta file to split

optional arguments:
  -h, --help            show this help message and exit
  -p PARTS, --parts PARTS
                        Number of parts to slice the file [10]
```
---


## Multifasta general stats 

If you have a directory contanining fasta files (fa: either faa or fna) compute several stats, that are important when describing MAGs [See Table 1 Preprint De Anda et al., 2020](https://www.researchsquare.com/article/rs-39998/v1)


```bash

for i in *.fa; do seqkit stat $i  >> stats; done
for i in *.fa ; do perl gc.pl $i >$i.gc.tab ; done

#Sum the scaffold GC and get the average 

for i in *.tab; do awk '{sum+= $2; n++ } END { if (n > 0) print sum / n; }' $i > $i.GC.average ; done
 
```

---
## Generate sequence length file from multifasta

**Option 1 awk**

Obtained from [here](https://www.danielecook.com/generate-fasta-sequence-lengths/)

```bash
cat file.fa | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' 
```
```bash
awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' file.fa
```
**Option 2 [Seqkit](https://bioinf.shenwei.me/seqkit/)**

```bash
seqkit fx2tab --length --name --header-line file.fa >> file.lenght 
```

**Option3 [samtools](http://www.htslib.org/)**

```bash
samtools faidx file.fa  |  cut -f1-2 file.fa.fai > file.lenght.tab
```


## Convert fasta into 1 lners 

```
From this 

> header 1
ATGCAATGCATG
ATGCCCGGTAGT
TTATAGAGATAG

to this 

> header 1 
ATGCAATGCATGATGCCCGGTAGTTTATAGAGATAG
```

```perl
perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (sort(keys(%fa))){ print "$s\n$fa{$s}\n" }}' file.fa > file1ne.fa 
```

---
## Average length of multifasta

```perl
perl -lne 'if(/^(>.*)/){$h=$1}else{$fa{$h}.=$_} END{ foreach $h (keys(%fa)){$m+=length($fa{$h})}; printf("%1.0f\t",$m/scalar(keys(%fa))) }' file.fa
```
---


## Keep sequences of certain lenght 

In this case we are keeping sequences >100 bp

```perl
perl -lne 'if(/^(>.*)/){ $head=$1 } else { $fa{$head} .= $_ } END{ foreach $s (keys(%fa)){ print "$s\n$fa{$s}\n" if(length($fa{$s})>100) }}' file.fa > file100.fa
```
---

## Change the headers with the name of the name of fasta file

```perl
 perl -lne 'if(/^>(\S+)/){ print ">$ARGV $1"} else{ print }' file.fa > file_renamed.fa
```

---

## Histogram of total number of sequences in a large genomic dataset 

<img src="https://valdeanda.github.io/Useful_scripts/seq.png" width="400" height="300" align="right">

Let's suppose that you have thousands of genomes and you want to compare the total number of sequences in your genomic dataset. If all your genomes are  either .faa or .fna extension, you can use the following one-line command to count the number of sequences and generate a  histogram.  You can change the figure to pdf, just change pdf("seq.pdf");


```bash
grep -c ">" *.faa  | sed 's/:/\t/g' | cut -f 2 | Rscript -e 'data=abs(scan(file="stdin")); png("seq.png"); hist(data,xlab="secuences")'
```

---


## Remove fasta sequences from list of headers

It requires a list of headeres to remove from a fasta file  

```bash
python3 remove_sequences.py file.fa sequence_to_remove.txt > file_filtered.fa 
```
---

## Extract fasta sequences from list of headers

**Option 1 [pullseq](https://github.com/bcthomas/pullseq)**

```bash
pullseq -i file.fa -n  sequences_to_extract.txt > extracted_sequences.fa
```

**Option 2 [samtools](http://www.htslib.org/)**

```bash
cat sequences_to_extract.txt  | xargs -n 1 samtools faidx file.fa >> extracted_sequences.fa 
```
---
## Description of PFAM identifiers

1. Create a list of PFAM identifiers 

```bash
head identifiers.txt

PF13243
PF13249
PF02458
```

2. Run the following commands, originally created by [Dr. Carlos Cantalapiedra](https://github.com/Cantalapiedra) and incorporated in [MEBS](https://github.com/valdeanda/mebs)

```bash
cat identifieres.txt | while read  pfam; do
desc=$(curl http://pfam.xfam.org/family/"$pfam"/desc | head -1);
printf "$pfam\t";
printf "$desc\n";
done 2> /dev/null \
> identifiers.desc.tab
```


## Join files 

Take a column of 1 file and another column from another file and create a new file with those columns
No need for matching column 

```bash
paste <(awk '{print $1}' file1.txt ) <(awk '{print $2}' file2.txt ) > file3.txt
```

## Download genomes from a list ("ftp_GCA_download.txt") of ftp links

```bash
for next in $(cat ftp_GCA_download.txt); do wget  "$next"; done

or you can do wget -i ftp_GCA_download.txt
```
## Print duplicates in a column ($1 = column 1)

```bash
awk 'x[$1]++ == 1 { print $1 " is duplicated"}'
```

## Cut columns from 1 file and create a new file with those columns

```bash
cut -f2,3 file1.txt > file2.txt
```

## Search "pattern" and add "replace_pattern" at the beginning of the line (^)

```bash
sed '/pattern/ s/^/replace_pattern/' file.txt
```

## Search "pattern" and replace the 2nd occurrence of it (/2') with "replace_pattern"

```bash
sed 's/pattern/replace_pattern/2' file.txt
```

## Verifying empty columns 

In a file of 2 columns, if 2nd column of file is blank, print 1st column followed by "Your Words", otherwise print 1st and 2nd column, create new file of all this output

```bash
awk '{if (!$2) {print $1,"YourWords"} else {print $1, $2}}'  > file.tsv
```

## Download genomes from ncbi 
 

**Genome browse overview***
https://www.ncbi.nlm.nih.gov/genome/browse/#!/overview/

**Genbank assembly summary file

```bash
wget htp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt
```
**Get the complete and latest genomes from assembly summary genbank**

```bash
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary_genbank.txt
```


