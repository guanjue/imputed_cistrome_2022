args = commandArgs(trailingOnly=TRUE)

input_file = args[1]
output_file = args[2]

d = read.table(input_file, header=F, sep='\t')

dfdr = round(p.adjust(10^(-d[,4]), 'fdr'), 5)

write.table(cbind(d[,1:3], round(d[,4],5)/16)[dfdr<=0.01,], output_file, sep='\t', quote=F, col.names=F, row.names=F)
