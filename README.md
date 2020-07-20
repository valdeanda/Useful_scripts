# Useful_scripts

This repository contains scripts that I developed/modified from  my collaborator [Cesar Poot](https://github.com/acpooth) or useful tips/tricks I have picked up.


## Bubble plot

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
 
 Running Bubble plot with example data
 
 ```bash
  python3 bubble_chart.py  data_bubbleplot.tab -f  png -r 300
 ```

<img src="https://valdeanda.github.io/Useful_scripts/data_bubbleplot.tab_bubbleplot.png" width="400" height="300" align="right">
