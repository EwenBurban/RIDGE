args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
gff_file = args['gff_file']
contig_file = args['contig_file']
output = args['output']
window_size = as.numeric(args['window_size'])

gff = read.table(gff_file,h=F,sep='\t',quote="")[,1:5]
colnames(gff) = c('chr','source','type','start','end')
gff = subset(gff,gff$type=='exon')
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

contig_res_list = list()
for (i in 1:nrow(contig_data)){
	contig_length = contig_data[i,'contig_length']
	contig_name = contig_data[i,'contig_name_gff']
	contig_name_true = contig_data[i,'contig_name']
	contig_gff = subset(gff,subset=gff$chr == contig_name)
	win_seq = dna_window(1,contig_length,size=window_size)	
	win_res = lapply(win_seq,function(x,...){
			sub_gff = subset(contig_gff,subset=contig_gff['end'] > x[,'start'] & contig_gff['start'] < x[,'end'])
			if(nrow(sub_gff)==0){
				coding_rate=0
			}else{
				# correct for exon that could be paritaly out of the window
				if(any(sub_gff[,'end'] > x[,'end'])){sub_gff[which(sub_gff[,'end'] > x[,'end']),'end'] = x[,'end'] }
				if(any(sub_gff[,'start'] < x[,'start'])){sub_gff[which(sub_gff[,'start'] < x[,'start']),'start'] = x[,'start'] }

				# check if exon are correctly ordered
				diff = (sub_gff['end'] - sub_gff['start'])
				if (any(diff < 0)){print('your gff file is not correctly ordened (start pos that are after end pos)');q()}

				# correct for exon that could be partialy recouvrent with others
				sub_gff = sub_gff[order(sub_gff[,'start']),]
				v = as.vector(t(sub_gff[,c('start','end')]))
				dist =  v[2:length(v)] - v[1:(length(v)-1)]
				while(any(dist<0)){
					xkill = which(dist<0)[1]/2 +1 
					if(sub_gff[xkill,'end'] > sub_gff[xkill-1,'end']){sub_gff[xkill -1,'end'] = sub_gff[xkill,'end']}
					sub_gff = sub_gff[-xkill,]
					sub_gff = sub_gff[order(sub_gff[,'start']),]
					v = as.vector(t(sub_gff[,c('start','end')]))
					dist =  v[2:length(v)] - v[1:(length(v)-1)]
				}
				diff = (sub_gff['end'] - sub_gff['start'])
				coding_rate = sum(diff)/(x[,'end'] - x[,'start'])
			}
			res = c('chr'=contig_name_true,'start'=x[,'start'],'end'=x[,'end'],'conding_rate'=coding_rate)
			return(res)
		})
	contig_res_list[[i]] = do.call(rbind,win_res)
	print(contig_name)
}

full_res = do.call(rbind,contig_res_list)
colnames(full_res) = c('chr','start','end','coding_rate')
write.table(full_res,file=output,row.names=F,sep='\t')

