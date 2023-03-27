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
if(homo_rec!='True'){
	rho_file = args['rho_map']
	rho_map = read.table(rho_file,h=T,sep='\t')}

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

rec_rate_list = list()
for (i in 1:nrow(contig_data)){
	contig_name = contig_data[i,'contig_name']
	index = contig_data[i,'index']
	contig_length = contig_data[i,'contig_length']
	win_seq = dna_window(1,contig_length,size=window_size)
	y = 0
	tmp = lapply(win_seq,function(x,...){
		mean_point = mean(x)
		if(homo_rec=='True'){rec_rate=homo_rec_rate} else {
			rec_rate = subset(rho_map,subset=chr == index & start <= mean_point & end >= mean_point)['r'][[1]]
			if(identical(rec_rate,numeric(0))){rec_rate=NA}
			if(length(rec_rate)>1){rec_rate=mean(rec_rate,na.rm=T)}}
		res = c('chr' = contig_name, 'start' = x[,'start'] , 'end' = x[,'end'],'rec.rate' = rec_rate)
		return(res)
		})
	rec_rate_list[[i]] =  do.call(rbind,tmp)
	print(contig_name)
}
rec_rate_full = do.call(rbind,rec_rate_list)
colnames(rec_rate_full) = c('chr','start','end','rec_rate')
write.table(rec_rate_full,file=output,sep='\t',row.names=F)
