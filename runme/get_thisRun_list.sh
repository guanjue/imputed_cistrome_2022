target_list=$1
test_ct=$2
reliable_list=$3

### extracting training and testing cell types
rm -f TF_train_test.thisRun.list.txt
for tf_i in $(cat $target_list)
do
echo $tf_i
	cat $reliable_list | awk -F '\t' -v OFS='\t' -v tf_i=$tf_i -v test_ct=$test_ct '{if ($1==tf_i) print $0, toupper(test_ct)}' >> TF_train_test.thisRun.list.txt
done


