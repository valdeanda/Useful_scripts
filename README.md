# Useful_scripts

This repository contains scripts that I've developed with my collaborator [Cesar Poot](https://github.com/acpooth) and useful tips/tricks/scripts I have picked up over the years.

# Dependencies

Before using any of the scripts make sure that you have install the following dependencies

```bash
sudo apt-get install python3 python3-pip python3-matplotlib \
ipython3-notebook python3-mpltoolkits.basemap
sudo pip3 install -U pip
sudo -H pip3 install --upgrade pandas numpy scipy seaborn 
sudo -H pip3 install -U scikit-learn

```

---

## Bubble plot

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

# Split fasta 

Requires biopython

`` bash
pip3 install biopython
```

Script that is useful if you have a large fasta file and you want to split it into small files of the same size 

```bash
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


