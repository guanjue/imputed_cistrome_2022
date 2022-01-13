# imputed_cistrome
## about
imputed_cistrome predicts transcription factors binding sites based on chromatin accessibility and DNA sequence motif features. This project is still in progress.
## Dependences
* snakemake
* numpy
* pandas
* h5py
* pytorch
* pybedtools
* fire
* BeautifulSoup
* pybedtools
* pyBigWig
* scikit-learn
* giggle
## Usage
### modify the config file based on your environment
There are 2 example yaml config files (for human and mouse respectively). Please use them as templates and modify accordingly. After modifying the config file, run
```shell
python yaml_config_cistrome.py
```
### modify 'Snakefile' based on your target output
Here are some examples:
#### 1. train the models for a list of TFs
```python
list_tf_names = ['CTCF', 'ESR1', 'AR']

CELL_LINE, TF_NAME, CHROM_SET_NAME = tf_name_to_train(list_tf_names)

rule all:
    input:
        expand(
             "train/{tf_name}/models/{tf_name}.{training_cell_line}.{training_chrom_set_name}_model.pkl",
             zip,
             tf_name=TF_NAME,
             training_cell_line=CELL_LINE,
             training_chrom_set_name=CHROM_SET_NAME
         )
```
#### 2. train a TF model and predict the TF binding sites for some other cell type
```python
# the TF to be predicted
tf_name = 'CTCF'
# the cell type to be predicted, use one of the cell types defined in 'cell_types' variable in the yaml config file
cell_type = '48'

rule all:
    input:
        expand(
              "train/NR3C1/predictions/{tf_name}.{cell_type}.{chrom}_preds.h5",
              tf_name = tf_name,
              cell_type = cell_type,
              chrom = config['chrom_all']
        )
```
### execute the pipeline
It is recommended to run the script on a cluster.
First, dry run the pipeline to check whether there are problems:
```shell
snakemake -n
```
Then, submit jobs by:
```shell
snakemake -j JOB_NUM  --cluster-config cluster.yaml --latency-wait 60 --cluster "sbatch -J {cluster.job_name} -p {cluster.partition} -n {cluster.n}  -t {cluster.time} --mem {cluster.mem_mb}"
```
