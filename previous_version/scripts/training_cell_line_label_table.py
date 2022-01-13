import os
import yaml
import pandas as pd

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper

with open("../config.yaml", "r") as infile:
    config = yaml.load(infile, Loader=Loader)

# TF=[]
# for filename in os.listdir("./"):
#     if filename.endswith(".gz"):
#         TF.append(filename.split(".")[0])

data_dir = '/n/scratchlfs/xiaoleliu_lab/Jingyu/impute_cistrome/ENCODE_DREAM'


os.makedirs(config['training_cell_types_regions_label_path'], exist_ok=True)
os.makedirs(config['test_cell_types_regions_label_path'], exist_ok=True)

for tf_name in config['tf_names']:
    train_path_file = "ChIPseq/labels/%s.train.labels.tsv.gz" % tf_name
    test_heldout_chr_path_file = "ChIPseq/heldout_chr/training_celltypes/labels/%s.train.labels.tsv.gz" % tf_name
    df_train = pd.read_csv("%s/%s" % (data_dir, train_path_file), sep="\t", header=0)
    df_test_heldout_chr = pd.read_csv("%s/%s" % (data_dir, test_heldout_chr_path_file), sep="\t", header=0)
    df_train_all = pd.concat([df_train, df_test_heldout_chr], ignore_index=True)
    df_train_all.to_csv("%s/%s.%s" % (
        config['training_cell_types_regions_label_path'], tf_name, config['training_cell_types_regions_label_name']),
                        sep="\t", header=True, index=False)
    print(tf_name)


os.system('cp -r %s/%s %s' % (data_dir, "eva/labels/final/*", config['test_cell_types_regions_label_path']))


import os

for filename in os.listdir("./"):
    os.renames(filename, "%s_all_cell_types.h5" % filename.split(".")[0])



with open("config.yaml", "r") as infile:
    config = yaml.load(infile, Loader=Loader)

import os
import yaml
import pandas as pd

try:
    from yaml import CLoader as Loader, CDumper as Dumper
except ImportError:
    from yaml import Loader, Dumper


with open("config.yaml", "r") as infile:
    config = yaml.load(infile, Loader=Loader)

# path="label/test/"
path="label/train/"
for filename in os.listdir(path):
    df=df_all_regions_label = pd.read_csv(
        "%s/%s" % (path, filename),
        sep="\t", header=0, nrows=10)
    for cell_line in df.columns[3:]:
        if cell_line not in config["cell_types"]:
            print(filename,cell_line)

('CTCF.train.labels.tsv', 'IMR-90')
('MAFK.train.labels.tsv', 'IMR-90')
('CEBPB.train.labels.tsv', 'IMR-90')