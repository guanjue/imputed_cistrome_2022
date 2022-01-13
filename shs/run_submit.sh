#while read tf ct; do echo $tf; bash submit_getpks.sh $tf GM12878 K562 F 9; done < GM12878.txt

#bash submit.sh TAL1 K562_Erythroblast_Bone_Marrow None_Macrophage_Blood T

#cat /homes1/gxiang/projects/impute_cistrome/get_motif_difscore/reliable_TFs.txt | awk -F '\t' -v OFS='\t' -v test_ct=$test_ct '{print $1,$2,test_ct}' > TF_train_test.$test_ct.list.txt

time bash get_thisRun_list.sh

while read tf train_ct test_ct
do
	echo $tf
	bash submit.sh $tf $train_ct $test_ct T
done < TF_train_test.thisRun.list.txt
