#!/usr/bin/env Rscript

# Goal of the script
# This script aim to produce a file defining the genomic properties of each locus defined in {bedfile}
# This information are mandatory to simulate dataset containing locus with the same properties than observed dataset.

### gather command args and define variables
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
		y = vector()
		y[tmp[1]] = tmp[2]
		return(y)},USE.NAMES=F)
abc_data=read.table(args['abc_data'],h=T)
print('test')
Nref = as.numeric(args['Nref'])
window_size = as.numeric(args['window_size'])
popfile = read.table(args['popfile'],sep=',',h=T)
contig_data=read.table(args['contig_data'],h=T)
nameA = args['nameA']
nameB = args['nameB']
popA = popfile[,nameA]
popB = popfile[,nameB]
ploidy= as.numeric(args['ploidy'])
size_popA = length(popA[!is.na(popA)]) * ploidy
size_popB = length(popB[!is.na(popB)]) * ploidy
output = args['output']
nLoci=as.numeric(args['nLoci'])
output_bed=args['output_bed']
mu=as.numeric(args['mu'])
hetero_theta=args['hetero_theta']
print('full_load')
data=abc_data[sample(1:nrow(abc_data),size=nLoci),]
#generate  vectors of the same size than data, assuming that number of sample available is the same at each loci 
data$size_popA  = rep(size_popA,nrow(data))# generate a vector of the same size than data, assuming that number of sample available is the same at each loci 
data$size_popB  = rep(size_popB,nrow(data))
data$totpopsize = data$size_popA + data$size_popB
#### compute expected theta
if(hetero_theta=='True'){
data$theta=data$bialsite_avg / sum(1/(1:data$totpopsize[1]))
} else {
	data$theta=4 * Nref * mu * window_size
}
#### compute the rho value 
homo_rec=args['homo_rec']
homo_rec_rate=as.numeric(args['homo_rec_rate'])
if(homo_rec!='True'){ # If homo_rec==False, then the recombination rate is computed from a recombination map,
					  # Otherwise, if homo_rec==True, homo_rec_rate is applied homogeneously across all the genome.
	rho_file = args['rho_map']
	print('load rho map')
	rho_map = read.table(rho_file,h=T,sep='\t')
}


if(homo_rec=='True'){data$rho= homo_rec_rate * 4 * Nref * window_size} else {
	rec_vec=sapply(data$dataset,function(x,...){
					   str=unlist(strsplit(x,split=':|-'))
					   mean_point=mean(as.numeric(str[2:3]))
					   index=contig_data[which(contig_data$contig_name==str[1]),'index']
					   rec_rate = subset(rho_map,subset=chr == index & start <= mean_point & end >= mean_point)['r'][[1]] # gather the recombination rate inside the window
  					  if(identical(rec_rate,numeric(0))){rec_rate=0} # if recombination rate is equal to NA / miss (which is commonn at the start of chromosome, but depend of the recombination map), then replace by r=0
  					  if(length(rec_rate)>1){rec_rate=mean(rec_rate,na.rm=T)} # In some cases, the size of the window is higher than the window size of the recombination rate, so the recombination rate of the window is equal to the mean recombination rate inside the window
		return(rec_rate)
		})
	data$rho=rec_vec * 4 * Nref * window_size
	if(any(data$rho<0)){data$rho[which(data$rho<0)]=0}
}


data$locus_length = rep(window_size,nrow(data)) 
print(head(data))

# It keep only the necessary information for the rest of the pipeline and store it in {output} file
res = subset(data,select = c('locus_length','size_popA','size_popB','totpopsize','theta','rho'))
write.table(res,file=output,sep='\t',row.names=F)
bed=do.call(rbind,strsplit(data$dataset,split=':|-'))
colnames(bed)=c('chr','start','end')
write.table(bed,file=output_bed,sep='\t',row.names=F)




