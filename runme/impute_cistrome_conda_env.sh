# conda-4.11.0 update at 01/13/2022
conda update -n base conda

# conda install python3 and packages
conda create -n impute_cistrome_2022 python=3.9 r=4.1 fire pandas numpy scipy h5py matplotlib \
seaborn hyperopt pybedtools pyBigWig lightgbm=3.1.1 pyyaml

# other used default packages
gc os sys pickle functools timeit collections operator json

### activate conda and install two packages that is not available in conda
conda activate impute_cistrome_2022
pip install sklearn
pip install torch


### test run
cd /liulab/gjxiang/projects/test_impute_cistrome
git clone https://github.com/guanjue/impute_cistrome_2022.git
cp impute_cistrome_2022/runme/*sh test_GATA1/

cd test_GATA1
bash run_submit.done.sh






'''
conda create -n impute_cistrome_2022 python=3.9 fire=0.4.0 pandas=1.2.3 numpy=1.20.1 scipy=1.4.1 h5py=2.10.0 matplotlib=3.4.3 \
h5py=2.10.0 matplotlib=3.4.3 seaborn=0.11.1 \
hyperopt=0.2.5 pandas=1.2.3 numpy=1.20.1 scipy=1.4.1 pybedtools=0.8.2 \
pyBigWig=0.3.18 lightgbm=3.1.1


json=2.0.9
sklearn=0.24.1
torch=1.8.0
yaml=5.4.1
giggle=0.1.1
'''

