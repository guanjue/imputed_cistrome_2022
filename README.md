# impute_cistrome_GX

### Install conda env
### conda-4.11.0 update at 01/13/2022
```
conda update -n base conda
```

### conda install python3 and packages
```
conda create -n impute_cistrome_2022 python=3.9 r=4.1 fire pandas numpy scipy h5py matplotlib \
seaborn hyperopt pybedtools pyBigWig lightgbm=3.1.1 pyyaml
```

### activate conda and install two packages that is not available in conda
```
conda activate impute_cistrome_2022
pip install sklearn
pip install torch
```


### Test run
```
### git clone scripts
cd /liulab/gjxiang/projects/test_impute_cistrome
git clone https://github.com/guanjue/impute_cistrome_2022.git
```

### Step1: Enter working directory & copy input parameter files
```
mkdir test_GATA1
cp /liulab/gjxiang/projects/test_impute_cistrome/impute_cistrome_2022/runme/* test_GATA1/
```

### Step2: Set up parameters in 'run_submit.done.sh' in the 'test_GATA1/' folder ###
```
### The parameters in 'run_submit.done.sh' are as follows:
# set impute TF; training cell type; testing cell type; User defined sample id for DNase data
tf=GATA1
train_ct=K562
test_ct=NONE_B_LYMPHOCYTE_BLOOD
user_train_DNase_list=F
user_test_DNase_list=F

### set working directory
working_dir0='/liulab/gjxiang/projects/test_impute_cistrome/test_GATA1'
### set script directory
script_dir='/liulab/gjxiang/projects/test_impute_cistrome/imputed_cistrome_2022'

### set input file directory
DNase_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/DNase'
TF_json_folder='/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/TFs'
pksumit_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_summits/'
pk_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome/human_TF_peaks/'
motif_h5_folder='/liulab/gjxiang/projects/impute_cistrome/cistrome_impute_results_human/hdf5s/motif'

### DNase-seq metadata
DNase_info=$script_dir'/config_file/cistrome.metadata.Homo_sapiens.DNase.withheader.txt'
### TF ChIP-seq metadata
TF_info=$script_dir'/config_file/cistromeDB_human_TF_ChIP_samples_meta_info_peaks_2000_motifs_enrichment_FRiP_0.01_UDHS_0.7.xls'
### hg38 200bp bins with 50bp sliding window
bins=$script_dir'/config_file/hg38.200_50slide.bins.bed'
######
```

### Step3: Run imputed_cistrome_2022
```
# enter working directory
cd test_GATA1
bash run_submit.done.sh
```

### Step2.0: Get parameters using the following scripts
```
### Step2.0.1: Select imputed TFs, write them into the 'TF_select.tmp.txt' file
target_TF_list='TF_select.tmp.txt'

### Step2.0.2: Select impute cell types
test_ct='NONE_B_LYMPHOCYTE_BLOOD'

### Step2.0.3: Select reliable TF list
reliable_TFs_list=/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt
time bash get_thisRun_list.sh $target_TF_list $test_ct $reliable_TFs_list

### Step2.0.4: Get impute_cistrome pipeline input parameters: impute TF; training cell type; testing cell type
cat TF_train_test.thisRun.list.txt

```
