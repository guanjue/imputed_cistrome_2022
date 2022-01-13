# impute_cistrome_GX

### setting up conda env
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


### test run
```
### git clone scripts
cd /liulab/gjxiang/projects/test_impute_cistrome
git clone https://github.com/guanjue/impute_cistrome_2022.git

### Step1: enter working directory & copy input parameter files
mkdir test_GATA1
cp /liulab/gjxiang/projects/test_impute_cistrome/impute_cistrome_2022/runme/* test_GATA1/

### Step2: Setting up parameters in 'run_submit.done.sh'
cd test_GATA1

### Step3: Run imputed_cistrome_2022
bash run_submit.done.sh
```

### Step2.0: Get parameters using the following scripts
```
### Step2.0.1: select imputed TFs, write them into the 'TF_select.tmp.txt' file
target_TF_list='TF_select.tmp.txt'

### Step2.0.2: select impute cell types
test_ct='NONE_B_LYMPHOCYTE_BLOOD'

### Step2.0.3: filter to select reliable TF list
reliable_TFs_list=/liulab/gjxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt
time bash get_thisRun_list.sh $target_TF_list $test_ct $reliable_TFs_list

### Step2.0.4: get impute_cistrome pipeline input parameters: impute TF; training cell type; testing cell type
cat TF_train_test.thisRun.list.txt

```
