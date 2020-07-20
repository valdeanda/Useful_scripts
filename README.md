# Useful_scripts

This repository contains scripts that I developed/modified from  my collaborator [Cesar Poot](https://github.com/acpooth) or useful tips/tricks I have picked up.

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
