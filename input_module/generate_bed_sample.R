#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
contig_file = args['contig_file']
window_size = as.numeric(args['window_size'])
nLoci_per_chr = as.numeric(args['nLoci_per_chr'])
output = args['output']
contig_data = read.table(contig_file,h=T,sep='\t')

dna_window = function(start,end,size,step=NULL){
	pos=start
	if(is.null(step)){step=size}
	a=vector()
	b=vector()
	while(pos + size -1 <= end){
		a[length(a) + 1] = pos
		b[length(b) + 1] = pos + size -1
		pos = pos + step
	}
	if(b[length(b)] < end){
		a[length(a) + 1] = b[length(b)] + 1
		b[length(b) + 1] = end
	}
	y=cbind(a,b)
	colnames(y) = c('start','end')
	y = split.data.frame(y,seq(nrow(y)))
	return(y)
}
sel_win_list = list()
for (i in 1:nrow(contig_data)){
	contig_name = contig_data[i,'contig_name']
	contig_length = contig_data[i,'contig_length']
	win_seq = dna_window(1,contig_length,size=window_size)
	if(nLoci_per_chr == -1){sel_win = win_seq} else{
	sel_win = sample(win_seq,nLoci_per_chr,replace=F)}
	sel_win = do.call(rbind,sel_win)
	print(nrow(sel_win))
	sel_win = cbind(rep(contig_name,nrow(sel_win)),sel_win)
	sel_win_list[[i]] = sel_win
}

bed = do.call(rbind,sel_win_list)
colnames(bed)= c('chr','start','end')
write.table(bed,file=output,sep='\t',row.names=F)

