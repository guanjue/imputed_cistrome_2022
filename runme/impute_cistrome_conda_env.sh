# conda-4.11.0 update at 01/13/2022
conda update -n base conda

# conda install python3 and packages
conda create -n impute_cistrome_2022 python=3.9 r=4.1 fire pandas numpy scipy h5py matplotlib \
seaborn hyperopt pybedtools pyBigWig lightgbm=3.1.1 pyyaml

### activate conda and install two packages that is not available in conda
conda activate impute_cistrome_2022
pip install sklearn
pip install torch


### test run
cd /liulab/gjxiang/projects/test_impute_cistrome
git clone https://github.com/guanjue/impute_cistrome_2022.git
cp impute_cistrome_GX/runme/*sh test_GATA1/

cd test_GATA1
bash run_submit.done.sh



gc os sys pickle functools timeit collections operator json


conda create -n impute_cistrome_2022 python=3.9 fire=0.4.0 pandas=1.2.3 numpy=1.20.1 scipy=1.4.1 h5py=2.10.0 matplotlib=3.4.3


h5py=2.10.0 matplotlib=3.4.3 seaborn=0.11.1 \
hyperopt=0.2.5 pandas=1.2.3 numpy=1.20.1 scipy=1.4.1 pybedtools=0.8.2 \
pyBigWig=0.3.18 lightgbm=3.1.1


conda activate impute_cistrome_2022


json=2.0.9
sklearn=0.24.1
torch=1.8.0
yaml=5.4.1
giggle=0.1.1


conda install -n impute_cistrome_2022 r-devtools



PackagesNotFoundError: The following packages are not available from current channels:

  - json=2.0.9
  - sklearn=0.24.1
  - torch=1.8.0
  - python=3.85
  - yaml=5.4.1
  - giggle=0.1.1

Current channels:

  - https://repo.anaconda.com/pkgs/main/linux-64
  - https://repo.anaconda.com/pkgs/main/noarch
  - https://repo.anaconda.com/pkgs/r/linux-64
  - https://repo.anaconda.com/pkgs/r/noarch
  - https://conda.anaconda.org/conda-forge/linux-64
  - https://conda.anaconda.org/conda-forge/noarch
  - https://conda.anaconda.org/bioconda/linux-64
  - https://conda.anaconda.org/bioconda/noarch
  - https://conda.anaconda.org/liulab-dfci/linux-64
  - https://conda.anaconda.org/liulab-dfci/noarch

To search for alternate channels that may provide the conda package you're



gc
os
sys
pickle
fire
h5py
matplotlib
seaborn
hyperopt
sklearn
yaml
pandas
numpy
scipy
functools
pybedtools
json
torch
timeit
lightgbm
collections
operator
pyBigWig
#assignment
giggle



import gc
import os
import pickle

import fire
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
from hyperopt.fmin import generate_trials_to_calculate
from sklearn.preprocessing import StandardScaler


import os
import yaml
import pandas as pd

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


import fire
import pandas as pd
import numpy as np
from scipy.stats import gamma
from scipy.stats import poisson

#!/usr/bin/env python3
import os

import h5py

from functools import partial
import pybedtools
import json
import numpy as np
import torch
from torch.multiprocessing import Pool, Process
import gc


from hyperopt import STATUS_OK
from hyperopt import hp
from timeit import default_timer as timer
import numpy as np
import lightgbm as lgb
from hyperopt import tpe
from hyperopt import Trials

from hyperopt import fmin

import collections
import warnings
from operator import gt, lt
import numpy as np
import sys
import lightgbm as lgb
from lightgbm.callback import EarlyStopException
from lightgbm.callback import _format_eval_result
from lightgbm.compat import range_

import os
import fire
import h5py
import pickle
import pyBigWig
import pandas as pd
import numpy as np
from sklearn.preprocessing import QuantileTransformer

from assignment import linear_assignment
from numpy import zeros,abs,exp,sort,zeros_like,argmin,delete,mod,mean,median,dot,array,sum,asarray,sign
from numpy.random import permutation
from operator import itemgetter








conda install -n impute_cistrome_2022 -c bioconda bioconductor-fabia

conda install -n impute_cistrome_2022 r-Rfast r-WGCNA r-dynamicTreeCut

conda install -n impute_cistrome_2022 -c bioconda bioconductor-GO.db

conda install -n impute_cistrome_2022 r-devtools

