train_ct=$1
test_ct=$2
DNase_info=$3
DNase_h5_folder=$4
TF=$5

echo $train_ct
cut -f1,7,8,9,27 $DNase_info | awk -F '\t' -v OFS='\t' -v celltype=$train_ct '{if ((toupper($2)=="NONE") && (toupper($2"_"$3"_"$4)==toupper(celltype)) && ($5>0.05) ) print $0; else if ((toupper($2)!="NONE") && (toupper($2)==toupper(celltype)) && ($5>0.05)) print $0}' > $train_ct'.id.list.txt1'
cut -f1,7,8,9,27 $DNase_info | awk -F '\t' -v OFS='\t' -v celltype=$test_ct '{if ((toupper($2)=="NONE") && (toupper($2"_"$3"_"$4)==toupper(celltype)) && ($5>0.05) ) print $0; else if ((toupper($2)!="NONE") && (toupper($2)==toupper(celltype)) && ($5>0.05)) print $0}' > $test_ct'.id.list.txt1'

### rm not exsit files
rm -f $train_ct'.id.list.txt'
while read id cl ct tis fdr
do
	h5file=$DNase_h5_folder'/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr1.'$id'.h5'
	if test -f "$h5file"; then
		printf "%s\t%s\t%s\t%s\t%s\n" $id $cl $ct $tis $fdr >> $train_ct'.id.list.txt'
	fi
done < $train_ct'.id.list.txt1'

rm -f $test_ct'.id.list.txt'
while read id cl ct tis fdr
do
	h5file=$DNase_h5_folder'/DNASE_bam_5_mer_variable_bp_all_samples_lightGBM_chr1.'$id'.h5'
	if test -f "$h5file"; then
		printf "%s\t%s\t%s\t%s\t%s\n" $id $cl $ct $tis $fdr >> $test_ct'.id.list.txt'
	fi
done < $test_ct'.id.list.txt1'

rm -f $train_ct'.id.list.txt1'
rm -f $test_ct'.id.list.txt1'

### filter top 10
train_num=$(wc -l $train_ct'.id.list.txt' |awk '{print $1}')
if (( $train_num > 10 ))
then
	sort -k5,5nr $train_ct'.id.list.txt' | head -10 | cut -f1 > $train_ct'.id.list.txt.tmp' && mv $train_ct'.id.list.txt.tmp' $train_ct'.id.list.txt'
else
	cut -f1 $train_ct'.id.list.txt' > $train_ct'.id.list.txt.tmp' && mv $train_ct'.id.list.txt.tmp' $train_ct'.id.list.txt'
fi



test_num=$(wc -l $test_ct'.id.list.txt' |awk '{print $1}')
if (( $test_num > 10 ))
then
        sort -k5,5nr $test_ct'.id.list.txt' | head -10 | cut -f1 > $test_ct'.id.list.txt.tmp' && mv $test_ct'.id.list.txt.tmp' $test_ct'.id.list.txt'
else
	cut -f1 $test_ct'.id.list.txt' > $test_ct'.id.list.txt.tmp' && mv $test_ct'.id.list.txt.tmp' $test_ct'.id.list.txt'
fi


