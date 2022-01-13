# impute_cistrome_GX

# setting up conda env

# conda-4.11.0 update at 01/13/2022
conda update -n base conda

# conda install python3 and packages
conda create -n impute_cistrome_2022 python=3.9 r=4.1 fire pandas numpy scipy h5py matplotlib \
seaborn hyperopt pybedtools pyBigWig lightgbm=3.1.1 pyyaml

### activate conda and install two packages that is not available in conda
conda activate impute_cistrome_2022
pip install sklearn
pip install torch



