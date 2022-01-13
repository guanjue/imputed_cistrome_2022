
reliable_list='/homes1/gxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt'
target_list='TF_select.txt'
test_ct='None_Macrophage_Blood'


rm -f TF_train_test.thisRun.list.txt
for tf_i in $(cat $target_list)
do
echo $tf_i
	cat $reliable_list | awk -F '\t' -v OFS='\t' -v tf_i=$tf_i -v test_ct=$test_ct '{if ($1==tf_i) print $0, toupper(test_ct)}' >> TF_train_test.thisRun.list.txt
done


