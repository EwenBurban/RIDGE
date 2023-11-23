#!/usr/bin/env Rscript

# Main goal of the script
# This script aim to create a file containing the recombination rate per base / per generation / per individual 
# at define windows. It can compute it considering an homogeneous recombination rate across the genome or
# considering an heterogeneous recombination rate computed from a user-provided recombination map.


### gather command args and define variables
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
contig_file = args['contig_file']
output = args['output']
window_size = as.numeric(args['window_size'])
contig_data = read.table(contig_file,h=T,sep='\t')
homo_rec=args['homo_rec']
homo_rec_rate=as.numeric(args['homo_rec_rate'])
if(homo_rec!='True'){ # If homo_rec==False, then the recombination rate is computed from a recombination map,
					  # Otherwise, if homo_rec==True, homo_rec_rate is applied homogeneously across all the genome.
	rho_file = args['rho_map']
	rho_map = read.table(rho_file,h=T,sep='\t')}

# This function is identical at the dna_window function from generate_bed_sample.R
#This function aim to generate  windows of {window_size} bp between two point (start and end)
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

# The following process, is applied on each chromosomes/contig separatly. 
#It generate genomic window in the same manner than in generate_bed_sample.R, in order to be merged with generate_bed_sample output file.
# When all window are generated, recombination rate is computed. 
#									If homo_rec == True, the recombination rate for each window is equal to {homo_rec_rate}
#									If homo_rec == False, the recombination rate is the average recombination rate in the window.
											
rec_rate_list = list()
for (i in 1:nrow(contig_data)){
	contig_name = contig_data[i,'contig_name']
	index = contig_data[i,'index']
	contig_length = contig_data[i,'contig_length']
	win_seq = dna_window(1,contig_length,size=window_size) # generate window list
	y = 0
	tmp = lapply(win_seq,function(x,...){ # apply for each window
		mean_point = mean(x) # determine the center position of the window
		if(homo_rec=='True'){rec_rate=homo_rec_rate} else { 
			rec_rate = subset(rho_map,subset=chr == index & start <= mean_point & end >= mean_point)['r'][[1]] # gather the recombination rate inside the window
			if(identical(rec_rate,numeric(0))){rec_rate=NA} # if recombination rate is equal to NA / miss (which is commonn at the start of chromosome, but depend of the recombination map), then replace by r=0
			if(length(rec_rate)>1){rec_rate=mean(rec_rate,na.rm=T)}} # In some cases, the size of the window is higher than the window size of the recombination rate, so the recombination rate of the window is equal to the mean recombination rate inside the window
		res = c('chr' = contig_name, 'start' = x[,'start'] , 'end' = x[,'end'],'rec.rate' = rec_rate)
		return(res)
		})
	rec_rate_list[[i]] =  do.call(rbind,tmp) 
	print(contig_name)
}
rec_rate_full = do.call(rbind,rec_rate_list) # compile results from all chr/contig in one data.frame and write it in {output} file.
colnames(rec_rate_full) = c('chr','start','end','rec_rate')
write.table(rec_rate_full,file=output,sep='\t',row.names=F)
