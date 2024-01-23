require(bio3d)
require(stringr)
aa=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-")
level=factor(aa)
f=read.fasta("./data/gap_le_10.fasta")
allcount=table(factor(f$ali,levels = level))
p=allcount/sum(allcount)
len=length(f$ali[1,])
probs=c()
gap=c()
ps=c()

for (i in seq(1,len)){
  ps=rbind(ps,p)
  tbl=table(factor(f$ali[,i],levels = level))  
  prob=tbl/sum(tbl)
  probs=rbind(probs,prob)
  if(max(prob)==prob[21]){
    gap=c(gap,T)
  }else{
    gap=c(gap,F)
  }
}

q=as.matrix(probs[,]+0.00000001)
S2=apply(-q*log(q),1,sum)
resi=seq(1,len)
pkt=which(is.element(resi,seq(61,100)))
resi_pkt = resi[pkt]
res1=resi_pkt[order(S2[pkt])]
print(paste0(res1[1:10],collapse = "+"))
print(sort(res1[1:10]))
pkt=which(is.element(resi,seq(242,282)))
resi_pkt = resi[pkt]
res2=resi_pkt[order(S2[pkt])]
print(paste0(res2[1:10],collapse = "+"))
print(sort(res2[1:10]))
