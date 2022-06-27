#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
            y = vector()
            y[tmp[1]] = tmp[2]
            return(y)},USE.NAMES=F)
vcf=args['vcf']
output = args['output']
window_size  = as.numeric(args['window_size'])
popfile = args['popfile']
contig_file = args['contig_file']
nameA = args['nameA']
nameB = args['nameB']
safe_mode='OFF'
rollwindow=10000
## get  the number of samples
header=grep('#CHROM',readLines(vcf),value=T)
header=sub('\n','',header)
header=unlist(strsplit(header,split='\t'))
popfile = read.table(popfile,sep=',',h=T)
header=header[-match(setdiff(header,c(popfile[,nameA],popfile[,nameB])),header)]
popA = match(popfile[,nameA],header)
if(any(is.na(popA))){popA = popA[-which(is.na(popA))]}
popB = match(popfile[,nameB],header)
if(any(is.na(popB))){popB = popB[-which(is.na(popB))]}
name_pop = rep(NA,sum(length(popA) + length(popB)))
name_pop[popA] = nameA
name_pop[popB] = nameB
## get the number of pertinent K
library(pcadapt)
data = read.pcadapt(vcf,type='vcf')
pca  = pcadapt(data,K=20)
pvar = pca$singular.values^2
var_of_interest = which(pvar >=0.05) # Keep only the PC that explain more than 5% of the variance
var2keep = vector()
for( i in var_of_interest){
	score = pca$scores[,i]
	a = score[popA]
	b = score[popB]
	# Test if the PC is able to separate the two groups
	testA =prod(max(a,na.rm=T) - quantile(b,c(0,1),na.rm=T))
	testB =prod(max(b,na.rm=T) - quantile(a,c(0,1),na.rm=T))
	if(( testA >=0) & ( testB>=0 )){
		var2keep[length(var2keep)+1] = i
	}
}
if (length(var2keep)==0 & safe_mode=='ON'){
	print("Warning, no PC has been able to distinguish the two population")
	print("Maybe there is an error in the attribution of sample to population")
	pdf('control_pcadapt.pdf')
	plot(pca,option='scores',pop=name_pop)
	dev.off()
	q() 
}

if (length(var2keep)==0 & safe_mode=='OFF'){
	print("Warning, no PC has been able to distinguish the two population")
	print("Maybe there is an error in the attribution of sample to population")
	print("Safe_mode OFF, process continue using only the first axis, results migth be weirds")
	var2keep=1
}
## get the pvalues and link it with the locus position
pca = pcadapt(data,K=max(var2keep))
pvalues = pca$pvalues
vcf = read.table(vcf,sep='\t')[,c(1:2)]
colnames(vcf) = c('chr','pos')
full_data = cbind(vcf,pvalues)
rm('vcf')
## rollwindowing to get max local pvalues
full_data$log = -log10(full_data$pvalues)
 
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
contig_data = read.table(contig_file,h=T)

#scan all the genome to find local maximum of -log10(pvalues)
max_pos_chr = list()
for (i in 1:nrow(contig_data)){
	contig = contig_data[i,'contig_name']
	contig_length = contig_data[i,'contig_length']
	tmp_data = subset(full_data,subset = chr == contig)

	win_seq = dna_window(1,contig_length,size=rollwindow,step=rollwindow/2)
	max_pos = lapply(win_seq,function(x,...){
			y=subset(tmp_data,subset=tmp_data$pos >= x[,'start'] & tmp_data$pos <= x[,'end'])
			maxrow = which(y$log == max(y$log,na.rm=T))
			return(y[maxrow,])
			})
	max_pos = do.call(rbind,max_pos)
	max_pos = max_pos[match(unique(max_pos$pos),max_pos$pos),] # remove doublons
	max_pos_chr[[i]] = max_pos	
}
full_max_pos = do.call(rbind,max_pos_chr)
#get the windowed mean of -log10(pvalues) and pvalues
windowed_data_list = list()
for (i in 1:nrow(contig_data)){
	contig = contig_data[i,'contig_name']
	contig_length = contig_data[i,'contig_length']
	tmp_data = subset(full_max_pos,subset = chr == contig)
	win_seq = dna_window(1,contig_length,size=window_size)
	mean_pcadapt = lapply(win_seq,function(x,...){
			y=subset(tmp_data,subset=tmp_data$pos >= x[,'start'] & tmp_data$pos <= x[,'end'])
			mean_pca = mean(y$pvalues,na.rm=T)
			mean_log = mean(y$log,na.rm=T)
			y = c('chr'=contig,'start'=x[,'start'],'end'=x[,'end'],'pvalues'=mean_pca,'log10_pval'=mean_log)
			return(y)
			})
	windowed_data_list[[i]] = do.call(rbind,mean_pcadapt)

}
windowed_data = do.call(rbind,windowed_data_list)
colnames(windowed_data) = c('chr','start','end','pvalues','log_pvalues')
if(any(windowed_data[,'pvalues'] ==  "NaN")){
windowed_data = windowed_data[-which(windowed_data[,'pvalues']=='NaN'),]}

write.table(windowed_data,file=output,row.names=F,sep='\t')
