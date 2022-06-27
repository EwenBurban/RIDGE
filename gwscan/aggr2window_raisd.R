#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
input_dir = args['input_dir']
output = args['output']
pop = args['pop']
contig_file = args['contig_file']
window_size  =as.numeric( args['window_size'])

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


contig_data = read.table(contig_file,sep='\t',h=T)
windowed_data_list = list()
for (i in 1:nrow(contig_data)){
	contig = contig_data[i,'contig_name']
	contig_length = contig_data[i,'contig_length']
	tmp_data = read.table(file.path(input_dir,paste0('RAiSD_Report.',pop,'.',contig)),h=F)
	colnames(tmp_data) = c('pos','mu')
	win_seq = dna_window(1,contig_length,size=window_size)
	win_mu = lapply(win_seq,function(x,...){
						y=subset(tmp_data,subset=tmp_data$pos >= x[,'start'] & tmp_data$pos <= x[,'end'])
						mu = mean(y$mu,na.rm=T)
						yy = c('chr'=contig,'start'=x[,'start'],'end'=x[,'end'],'mu'=mu)
						return(yy)})
	windowed_tmp_data = do.call(rbind,win_mu)
	print(head(windowed_tmp_data))
	windowed_data_list[[i]] = windowed_tmp_data
}
windowed_data = do.call(rbind,windowed_data_list)
colnames(windowed_data) = c('chr','start','end',paste0('mu_',pop))
if(any(windowed_data[,paste0('mu_',pop)] == "NaN")){
windowed_data= windowed_data[-which(windowed_data[,paste0('mu_',pop)]=="NaN"),]}
print(head(windowed_data))
write.table(windowed_data,file=output,row.names=F,sep='\t')

	


