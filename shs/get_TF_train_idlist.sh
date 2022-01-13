TF=$1
train_ct=$2
TF_info=$3
pksumit_folder=$4


cut -f1,3,9,10,11 $TF_info | awk -F '\t' -v OFS='\t' -v celltype=$train_ct '{if ((toupper($3)=="NONE") && (toupper($3"_"$4"_"$5)==toupper(celltype))) print $0; else if ((toupper($3)!="NONE") && (toupper($3)==toupper(celltype))) print $0}' > $train_ct'.'$TF'.txt'

cat $train_ct'.'$TF'.txt' | awk -F '\t' -v OFS='\t' -v TF=$TF '{if ($2==TF) print $1}' > $train_ct'.'$TF'.train.id.list.txt1'

rm -f $train_ct'.'$TF'.train.id.list.txt'
while read id
do
	submitfile=$pksumit_folder'/'$id'_sort_summits.bed'
	if test -f "$submitfile"; then
		printf "%s\n" $id >> $train_ct'.'$TF'.train.id.list.txt'
	fi
done < $train_ct'.'$TF'.train.id.list.txt1'
#rm $train_ct'.'$TF'.train.id.list.txt1'
#rm $train_ct'.'$TF'.txt'

TF_train_num=$(wc -l $train_ct'.'$TF'.train.id.list.txt' |awk '{print $1}')
if (( $TF_train_num == 1 ))
then
	cat $train_ct'.'$TF'.train.id.list.txt' > $train_ct'.'$TF'.train.id.list.txt.tmp'
	cat $train_ct'.'$TF'.train.id.list.txt.tmp' >> $train_ct'.'$TF'.train.id.list.txt'
	rm $train_ct'.'$TF'.train.id.list.txt.tmp'
fi

