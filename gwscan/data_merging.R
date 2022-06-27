#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=T)
args = sapply(args,function(x){tmp = unlist(strsplit(x,split='='))
           y = vector()
           y[tmp[1]] = tmp[2]
           return(y)},USE.NAMES=F)
abc_file=args['abc']
rec_rate_file = args['rec_rate']
coding_rate_file = args['coding_rate']
sweep_A=args['sweep_A']
sweep_B = args['sweep_B']
adapt = args['pcadapt']
outfile=args['output']
## custom func ##
get_pos = function(str){
    tmp=unlist(strsplit(str,split=':'))
    chr=tmp[1]
    tmp_bis = unlist(strsplit(tmp[2],split='-'))
    start=as.integer(tmp_bis[1])
    end=as.integer(tmp_bis[2])
    res = c('chr'=chr,'start'=start,'end'=end)
    return(res)
}
### abcstat
data = read.table(abc_file,h=T,sep='\t')
pos = do.call(rbind,lapply(data$dataset,get_pos))
data = cbind(pos,data[,-which(colnames(data)=='dataset')])
data[,-1] = apply(data[,-1],2,as.numeric)
print(head(data))
### coding rate
if(is.na(coding_rate_file)==F){
coding_rate = read.table(coding_rate_file,h=T,sep='\t')
data = merge(data,coding_rate,by=c('chr','start','end'),all=T)
}
print(head(data))
### rec rate
if(is.na(rec_rate_file)==F){
rec_rate = read.table(rec_rate_file,h=T,sep='\t')
data = merge(data,rec_rate,by=c('chr','start','end'),all=T)
}
print(head(data))
### sweeps
if(is.na(sweep_A)==F){
sweep_A = read.table(sweep_A,h=T,sep='\t')
data = merge(data,sweep_A,by=c('chr','start','end'),all=T)
sweep_B = read.table(sweep_B,h=T,sep='\t')
data = merge(data,sweep_B,by=c('chr','start','end'),all=T)
}
print(head(data))
### adapt pvalues
if(is.na(adapt)==F){
adapt  =read.table(adapt,h=T,sep='\t')
data = merge(data,adapt,by=c('chr','start','end'),all=T)
}
print(head(data))
write.table(x=data,file=outfile,row.names=F,sep='\t')
