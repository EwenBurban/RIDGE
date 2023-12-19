#!/usr/bin/env Rscript

# Goal of the script
# This script aim to produce a list of locus of defined size of {window_size}. 
# This list might be an exhaustive list or a sample of {nLoci_per_chr} window per chromosome/contig

### gather command args and define variables
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
contig_file = args['contig_file']
window_size = as.numeric(args['window_size'])
output = args['output']
contig_data = read.table(contig_file,h=T,sep='\t')


# This function aim to generate  windows of {window_size} bp between two point (start and end)
# by default, the function produce non-overlapping windows (step=NULL). If step < {window_size}, it produce overlapping windows.
# If step > {window_size}, it produce distant windows
dna_window = function(start,end,size,step=NULL){
	#setup variable
	pos=start
	if(is.null(step)){step=size}
	a=vector()
	b=vector()

	# control the length between start and end point. If length is < to window size, it keep the entierity of the regions in one window.
	if(pos + size -1 > end){print(paste0('little chr of ',end));return(data.frame('start'=1,'end'=end))}

	# generate bounds of window by iterative process
	while(pos + size -1 <= end){
		a[length(a) + 1] = pos
		b[length(b) + 1] = pos + size -1
		pos = pos + step
	}
	if(b[length(b)] < end){ # when the while loop condition is broken, it take the rest of the chromosome and put it in a window
		a[length(a) + 1] = b[length(b)] + 1
		b[length(b) + 1] = end
	}
	# generate a data.frame to store start and end bound of windows
	y=cbind(a,b)
	colnames(y) = c('start','end')
	y = split.data.frame(y,seq(nrow(y)))
	return(y)
}


# This section apply the dna_window function to each chromosomes/contig mentioned in contig_data file. 
# It uses the size of the chromosome/contig as the end of the regions to split in windows. By default the first bases is at position 1
# If {nLoci_per_chr} is different from -1, it subsample {nLoci_per_chr} window per chromosome/contig 
sel_win_list = list()
for (i in 1:nrow(contig_data)){
	contig_name = contig_data[i,'contig_name']
	contig_length = contig_data[i,'contig_length']
	win_seq = dna_window(1,contig_length,size=window_size)
	sel_win = do.call(rbind,win_seq)
	sel_win = cbind(rep(contig_name,nrow(sel_win)),sel_win)
	sel_win_list[[i]] = sel_win
}
## produce the final result and store it in {output} file.

bed = do.call(rbind,sel_win_list)
colnames(bed)= c('chr','start','end')
write.table(bed,file=output,sep='\t',row.names=F)

