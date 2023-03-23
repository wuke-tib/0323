# AAComics

Omics data analysis reveals the system-level constraint on cellular amino acid composition

## Software

To create a stand-alone environment named AAComics with Python 3 and all the required package versions (especially for cobrapy is also available), run the following code:

```shell
$ conda create -n AAComics python=3
```
```shell
$ conda activate AAComics
```
```shell
$ pip install ipykernel  
$ python -m ipykernel install --user --name AAComics --display-name "AAComics"  
$ pip install pandas
$ pip install Bio
$ pip install matplotlib
$ pip install seaborn
$ pip install cobra
```
  You can read more about using conda environments in the [Managing Environments](http://conda.pydata.org/docs/using/envs.html) section of the conda documentation. 

## Steps to reproduce the main analysis in the publication

Typical results can be reproduced by executing the Jupyter Python notebooks:

+ protein_composition.py

  ——Functions definition file used for following analysis：

  #####       analysis_workflow.ipynb

  #####       draw_figure_in_article.ipynb


## Run command line

### analysis_workflow.py tag:0-原始数据处理，1-数据都为1，random-原始数据[最小值，最大值], mean-1000个随机df取均值, new-新的酵母文件分析
```
python ./script/analysis_workflow.py -i MG1655 W3110 Bacillus Yeast_single_cell Corynebacterium_RNA_seq BW25113 DH1 -o ./analysis_result/random_data/ -tag random -mean F

python ./script/analysis_workflow.py -i MG1655 W3110 Bacillus Yeast_single_cell Corynebacterium_RNA_seq BW25113 DH1 -o ./analysis_result/initial_data/ -tag 0 -mean F

python ./script/analysis_workflow.py -i MG1655 W3110 Bacillus Yeast_single_cell Corynebacterium_RNA_seq BW25113 DH1 -o ./analysis_result/data_1/ -tag 1 -mean F

python ./script/analysis_workflow.py -i W3110 Bacillus Yeast_single_cell Corynebacterium_RNA_seq BW25113 DH1 -o ./analysis_result/mean_data/ -tag random -mean T
```

python ./script/analysis_workflow.py -i scerevisiae -o ./analysis_result/initial_data/ -tag new -mean F