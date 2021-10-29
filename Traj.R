#!/usr/bin/Rscript
library(bio3d)
library(ncdf4)
library(igraph)
library(ggplot2)
library(gplots)
library(reshape2)
#######set parameters################
#setwd(commandArgs(TRUE)[1])

path<-getwd()
print(path)
files<-unlist(strsplit(path,'/'))
pdbname<-files[length(files)]
dir.create('analysis')
setwd('./analysis')
Traj_folder<-'/data/sheng/MD/Antibody_MD/'
path_ANARCI<-'/usr/bin/ANARCI'
path_TMalign<-'/home/sheng/work/software/bin/TMalign'

#######read trajectory################
pdb<-read.pdb('../Temp_MD.pdb')

pdb$atom[pdb$atom$resid %in% 'HIE',]$resid<-'HIS'
seleH<-atom.select(pdb,elety='CA',chain='H')
seleL<-atom.select(pdb,elety='CA',chain='L')
ca.Hpdb<-trim(pdb,inds=seleH)
ca.Lpdb<-trim(pdb,inds=seleL)

seqH<-paste(c(as.character(pdbseq(pdb,inds=seleH))),sep='',collapse='')
seqL<-paste(c(as.character(pdbseq(pdb,inds=seleL))),sep='',collapse='')

H_numbering<-system(paste(path_ANARCI,' --scheme c --sequence ',seqH,sep=''),intern=T)
L_numbering<-system(paste(path_ANARCI,' --scheme c --sequence ',seqL,sep=''),intern=T)
H_numbering<-H_numbering[-c(1:7,length(H_numbering),grep('-',H_numbering))]
L_numbering<-L_numbering[-c(1:7,length(L_numbering),grep('-',L_numbering))]
seleposeH<-ca.Hpdb$atom$resno[length(H_numbering)]#determine last residue of variable domain
seleposeL<-ca.Lpdb$atom$resno[length(L_numbering)]
seleposeHst<-ca.Hpdb$atom$resno[1]#determine last residue of variable domain
seleposeLst<-ca.Lpdb$atom$resno[1]

cn<-colorRampPalette(c('blue','gray90','red'))
threads<-12
####################Reading trajectories#############
if(commandArgs(TRUE)[1] %in% 'NAMD'){
dcd<-read.dcd('../production.dcd')
}else{
dcd<-read.ncdf('../production.mdcrd')
}

##############RMSF and RMSD for heavy and light chain alone############

for(i in c('H','L')){
	print(i)
	ca.inds <- atom.select(pdb, elety="CA",chain=i)
	segpdb<-trim(pdb,inds=ca.inds)
	seleposs<-min(segpdb$atom$resno)
	selepose<-seleposeH
	if(i =='L'){selepose<-seleposeL}
	posrange<-seleposs:selepose
	ca.inds <- atom.select(pdb, elety="CA",chain=i,resno=posrange)
	segpdbca<-trim(pdb,inds=ca.inds)
	segpdb<-trim(pdb,inds=ca.inds)

	xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz,  mobile.inds=ca.inds$xyz, ncore=threads)

	####
	print("fitting coordinates\n")
	pc <- pca.xyz(xyz[,ca.inds$xyz])
	rd <- rmsd(pdb$xyz[ca.inds$xyz], xyz[,ca.inds$xyz]) 
	print("calculating RMSD within trajectory")
	rdpair <- rmsd(xyz[,ca.inds$xyz],fit=TRUE, ncore=threads) 
	rf <- rmsf(xyz[, ca.inds$xyz])

	write.table(rd, paste("rmsd_",i,"_Fv.txt",sep=''), sep="\t")
	write.table(rdpair, paste("rmsd_pair_",i,"_Fv.txt",sep=''), sep="\t")


	pdf(paste('RMSD_pair',i,'_Fv.pdf',sep=''),width=12,height=12)
	heatmap.2(rdpair, dendrogram='row', Colv=F,trace='none',colRow=cn(nrow(xyz)),colCol=cn(nrow(xyz)),cexCol=0.09)
	dev.off()


	pdf(paste('RMSF_',i,'_Fv.pdf',sep=''))
	plot(rf, ylab="RMSF", xlab='', typ="l",xaxt='n')
	posname<-paste(segpdb$atom[segpdb$atom$elety %in% 'CA',]$chain,segpdb$atom[segpdb$atom$elety %in% 'CA',]$resno,sep='_')
	axis(1,at=seq(1,length(rf)+10,10),labels=posname[seq(1,length(rf)+10,10)],las=2)
	dev.off()

	pdf(paste('PCA_',i,'_Fv.pdf',sep=''))
	pc <- pca.xyz(xyz[,ca.inds$xyz])
	plot(pc, col=cn(nrow(xyz)),mar=c(4,5,1,1) ) 
	dev.off()



	pdf(paste('RMSD_frame_',i,'_Fv.pdf',sep=''))
	plot(rd, typ="l", ylab="RMSD", xlab="Frame No.") 
	points(lowess(rd), typ="l", col="red", lty=2, lwd=2) 
	dev.off()

	pdf(paste('RMSD_distribution',i,'_Fv.pdf',sep=''))
	hist(rd, breaks=40, freq=FALSE, main="RMSD Histogram", xlab="RMSD") 
	lines(density(rd), col="gray", lwd=3) 
	dev.off()


	pdf(paste('PCA_positional_contribution_',i,'_Fv.pdf',sep=''))
	plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", ann=F,axes=F, typ="l",col='blue',ylim=c(0,1)) #sse=ss,
	points(pc$au[,2], typ="l", col="green",xaxt='n') 
	points(pc$au[,3], typ="l", col="red",xaxt='n') 
	axis(1,at=seq(seleposs,selepose+10,10),labels=posname[seq(seleposs,selepose+10,10)],las=2)
	axis(2,at=seq(0,1,0.1),labels=seq(0,1,0.1),las=2)
	dev.off()

	p1 <- mktrj.pca(pc, pc=1, b=pc$au[,1], file=paste("pc1",i,"_Fv.pdb",sep=''))
	p2 <- mktrj.pca(pc, pc=2, b=pc$au[,2], file=paste("Bio3d_pc2",i,"_Fv.pdb",sep=''))
	p3 <- mktrj.pca(pc, pc=3, b=pc$au[,3], file=paste("Bio3d_pc3",i,"_Fv.pdb",sep=''))


	pdf(paste('correlation_matrix',i,'_Fv.pdf',sep=''))
	#ca.inds.H <- atom.select(pdb, elety="CA", resno=1:125)
	#xyz.H <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds.H$xyz,  mobile.inds=ca.inds.H$xyz)
	#cij.H<-dccm(xyz.H[,ca.inds.H$xyz])
	cij<-dccm(xyz[,ca.inds$xyz])
	plot(cij,col.regions=bluered(10),at=seq(-1,1,0.2),resno=segpdb)
	dev.off()

	write.table(cij,paste('correlation_matrix_',i,'.txt',sep=''))

}
######################analysis for Fab############################
print("working on Fab ")
H_sele<-atom.select(pdb,elety="CA",chain='H', resno=c(seleposeHst:seleposeH))
L_sele<-atom.select(pdb,elety="CA",chain='L', resno=c(seleposeHst:seleposeL))
ca.inds<-combine.select(sel1=H_sele,sel2=L_sele, operator='OR')

HL_sele<-atom.select(pdb,elety="CA",chain=c('H','L'))
consca.inds<-combine.select(sel1=HL_sele,sel2=ca.inds, operator='not')
Hall_sele<-atom.select(pdb,elety="CA",chain='H')
Lall_sele<-atom.select(pdb,elety="CA",chain='L')
Hconst<-combine.select(sel1=Hall_sele,sel2=H_sele, operator='not')
Lconst<-combine.select(sel1=Lall_sele,sel2=L_sele, operator='not')

for(hl in c('HL_Fv','H_constant','L_constant','HL_constant')){
	if(hl %in% 'HL_constant'){ca.inds=consca.inds}else 
	if(hl %in% 'H_constant'){ca.inds=Hconst}else 
	if(hl %in% 'L_constant'){ca.inds=Lconst}

print("fitting coordinates")
segpdb<-trim(pdb,inds=ca.inds)
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz,  mobile.inds=ca.inds$xyz, ncore=threads)
pc <- pca.xyz(xyz[,ca.inds$xyz])
rd <- rmsd(pdb$xyz[ca.inds$xyz], xyz[,ca.inds$xyz]) 
rf <- rmsf(xyz[, ca.inds$xyz])
rdpair <- rmsd(xyz[,ca.inds$xyz],fit=TRUE, ncore=threads) 

posHL<-segpdb$atom[segpdb$atom$eleno %in% ca.inds$atom,]$resno
rfd<-as.data.frame(cbind(posHL,rf))
colnames(rfd)<-c('pos','rmsf')
write.table(rfd, paste("rmsf_",hl,".txt",sep=''), sep="\t")
write.table(rd, paste("rmsd_",hl,".txt",sep=''), sep="\t")
write.table(rdpair, paste("rmsd_pair_",hl,".txt",sep=''), sep="\t")

pdf(paste('RMSF_',hl,'.pdf',sep=''))
plot(rf, ylab="RMSF", xlab='', typ="l",xaxt='n')
posname<-paste(segpdb$atom[segpdb$atom$elety %in% 'CA',]$chain,segpdb$atom[segpdb$atom$elety %in% 'CA',]$resno,sep='_')
axis(1,at=seq(1,length(posHL),10),labels=posname[seq(1,length(posHL),10)],las=2)
dev.off()

pdf(paste('RMSD_pair_',hl,'.pdf',sep=''),width=12,height=12)
heatmap.2(rdpair, dendrogram='row', Colv=F,trace='none',colRow=cn(1000),colCol=cn(1000),cexCol=0.09)
dev.off()

pdf(paste('PCA_',hl,'.pdf',sep=''))
pc <- pca.xyz(xyz[,ca.inds$xyz])
plot(pc, col=cn(nrow(xyz)),mar=c(4,5,1,1) ) 
dev.off()

pdf(paste('PCA_positional_contribution_',hl,'.pdf',sep=''))
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", ann=F,axes=F, typ="l",col='blue',ylim=c(0,1)) #sse=ss,
points(pc$au[,2], typ="l", col="green",xaxt='n') 
points(pc$au[,3], typ="l", col="red",xaxt='n') 
axis(1,at=seq(1,length(posHL),10),labels=posname[seq(1,length(posHL),10)],las=2)
axis(2,at=seq(0,1,0.1),labels=seq(0,1,0.1),las=2)
dev.off()

pdf(paste('RMSD_frame_',hl,'.pdf',sep=''))
plot(rd, typ="l", ylab="RMSD", xlab="Frame No.") 
points(lowess(rd), typ="l", col="red", lty=2, lwd=2) 
dev.off()

cij<-dccm(xyz[,ca.inds$xyz])
net<-cna(cij,cutoff.cij=0.6)
#plot(net,segpdb)
pdf(paste('correlation_matrix_',hl,'.pdf',sep=''))
plot(cij,col.regions=bluered(10),at=seq(-1,1,0.2),resno=segpdb)
dev.off()
#pymol(cij,segpdb,exefile=path_pymol,type='script',omit=0.3)

node.betweenness <- betweenness(net$network,normalized = T)
pdf(paste('correlation_betweenness_',hl,'.pdf',sep=''),width=8,height=4)
plot(node.betweenness, xlab="Residue No", ylab="Centrality", type="h")
dev.off()
write.pdb(segpdb, b=normalize.vector(node.betweenness), file=paste("cna_betweenness_",hl,".pdb",sep=''))
}



###################VH-VL angles and distance###############################
Hcoreset<-c(35,12,38,36,83,19,94,37,11,47,39,93,46,45,68,69,71,70,17,72,92,84,91,90,20,21,85,25,24,86,89,88,87,22,23)
Lcoreset<-c(44,19,69,14,75,82,15,21,47,20,48,49,22,81,79,80,23,36,35,37,74,88,38,18,87,17,86,85,46,70,45,16,71,72,73)
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz,  mobile.inds=ca.inds$xyz)

seleH<-atom.select(pdb,elety='CA',chain='H')
seleL<-atom.select(pdb,elety='CA',chain='L')
ca.Hpdb<-trim(pdb,inds=seleH)
ca.Lpdb<-trim(pdb,inds=seleL)

seqH<-paste(c(as.character(pdbseq(pdb,inds=seleH))),sep='',collapse='')
seqL<-paste(c(as.character(pdbseq(pdb,inds=seleL))),sep='',collapse='')

H_numbering<-system(paste(path_ANARCI,' --scheme c --sequence ',seqH,sep=''),intern=T)
L_numbering<-system(paste(path_ANARCI,' --scheme c --sequence ',seqL,sep=''),intern=T)
H_numbering<-H_numbering[-c(1:7,length(H_numbering),grep('-',H_numbering))]
L_numbering<-L_numbering[-c(1:7,length(L_numbering),grep('-',L_numbering))]

H_match<-data.frame(cbind(ca.Hpdb$atom$resno[1:length(H_numbering)],gsub('^[A-Z ]+([0-9]+) [A-Z ]+$','\\1',H_numbering,perl=T)))
L_match<-data.frame(cbind(ca.Lpdb$atom$resno[1:length(L_numbering)],gsub('^[A-Z ]+([0-9]+) [A-Z ]+$','\\1',L_numbering,perl=T)))

H_corepos<-H_match[H_match$X2 %in% Hcoreset,]$X1
L_corepos<-L_match[L_match$X2 %in% Lcoreset,]$X1

matri<-matrix(nrow=nrow(dcd),ncol=6)
Hs<-atom.select(pdb,resno=as.numeric(as.vector(H_corepos)),elety='CA',chain='H')
Ls<-atom.select(pdb,resno=as.numeric(as.vector(L_corepos)),elety='CA',chain='L')

transformc<-function(matrix1,coord){
  x<-matrix1[1,1]+sum(matrix1[1,2:4]*coord)
  y<-matrix1[2,1]+sum(matrix1[2,2:4]*coord)
  z<-matrix1[3,1]+sum(matrix1[3,2:4]*coord)
  return(c(x,y,z));
}

normalize<-function(vector1){
  x<-sqrt(vector1[1]^2+vector1[2]^2+vector1[3]^2)
  return(vector1/x);
}

vcrossp <- function( a, b ) { 
  result=c(a[2] * b[3] - a[3] * b[2], a[3] * b[1] - a[1] * b[3], a[1] * b[2] - a[2] * b[1] )
  result 
} 

angle_dis<-function(Hpdb,Lpdb){
  system(paste(path_TMalign," ",Traj_folder,'/consensus_H.pdb ',Hpdb,' -m matrix_H.txt',sep=''))
  system(paste(path_TMalign," ",Traj_folder,'/consensus_L.pdb ',Lpdb,' -m matrix_L.txt',sep=''))
  
  pcH<-as.matrix(rbind(c(0.9525187, -0.1371821, 0.2718256),c(-0.117058, 0.659152, 0.7428432),c(-2.691829, -3.847092, 1.196887)))
  pcL<-as.matrix(rbind(c(-0.6193343, 0.639472, 0.4555223),c(0.5267385, 0.7686645, -0.362907),c(-3.702842, -0.6288583, -5.314558)))
  cH<-c(sum(pcH[,1]*c(-5,0.5,1)),sum(pcH[,2]*c(-5,0.5,1)),sum(pcH[,3]*c(-5,0.5,1)))
  cL<-c(sum(pcL[,1]*c(3,-1,1)),sum(pcL[,2]*c(3,-1,1)),sum(pcL[,3]*c(3,-1,1)))
  
  L1<-c(cL[1]+pcL[1,1],cL[2]+pcL[1,2],cL[3]+pcL[1,3])
  L2<-c(cL[1]+pcL[2,1],cL[2]+pcL[2,2],cL[3]+pcL[2,3])
  H1<-c(cH[1]+pcH[1,1],cH[2]+pcH[1,2],cH[3]+pcH[1,3])
  H2<-c(cH[1]+pcH[2,1],cH[2]+pcH[2,2],cH[3]+pcH[2,3])
  
  fileName <- "./matrix_H.txt"
  conn <- file(fileName,open="r")
  linn <-readLines(conn)
  transform_maxH<-matrix(c(as.numeric(unlist(strsplit(linn[3],'[\t ]+',perl=T))[3:6]),c(as.numeric(unlist(strsplit(linn[4],'[\t ]+',perl=T))[3:6]),c(as.numeric(unlist(strsplit(linn[5],'[\t ]+',perl=T))[3:6])))), nrow=3,byrow=T)
	close(conn)
  fileName <- "./matrix_L.txt"
  conn <- file(fileName,open="r")
  linn <-readLines(conn)
  transform_maxL<-matrix(c(as.numeric(unlist(strsplit(linn[3],'[\t ]+',perl=T))[3:6]),c(as.numeric(unlist(strsplit(linn[4],'[\t ]+',perl=T))[3:6]),c(as.numeric(unlist(strsplit(linn[5],'[\t ]+',perl=T))[3:6])))), nrow=3,byrow=T)
	close(conn)  
  cH<-transformc(transform_maxH,cH)
  cL<-transformc(transform_maxL,cL)
  L1<-transformc(transform_maxL,L1)
  L2<-transformc(transform_maxL,L2)
  H1<-transformc(transform_maxH,H1)
  H2<-transformc(transform_maxH,H2)
  
  C<-normalize(cH-cL)
  Cminus<--1*C
  L1<-normalize(L1-cL)
  L2<-normalize(L2-cL)
  H1<-normalize(H1-cH)
  H2<-normalize(H2-cH)
  dc<-sqrt(sum((cH-cL)^2))
  
  n_x=vcrossp(L1,C)
  n_y=vcrossp(C,n_x)
  tmpL=normalize(c(0,L1 %*% n_x,L1 %*%n_y))
  tmpH=normalize(c(0,H1 %*% n_x,H1 %*%n_y))
  
  HL<-acos(tmpL %*% tmpH)*(180/pi)
  if(HL>90){HL<-180-HL}
  if(HL>0){
    HL<--HL	
  }
  LC1<-acos(L1 %*% C)*(180/pi)
  LC2<-acos(L2 %*% C)*(180/pi)
  
  HC1<-acos(H1 %*% Cminus)*(180/pi)
  HC2<-acos(H2 %*% Cminus)*(180/pi)
  
  return(c(HL,HC1,HC2,dc,LC1,LC2))
}

for(i in 1:nrow(dcd)){
  write.pdb(pdb=trim(pdb,ind=Hs),file='H.pdb',xyz=xyz[i,Hs$xyz],print.segid=T)
  write.pdb(pdb=trim(pdb,ind=Ls),file='L.pdb',xyz=xyz[i,Ls$xyz],print.segid=T)
  P<-angle_dis('H.pdb','L.pdb')
  matri[i,]<-P
}
unlink(c('H.pdb', 'L.pdb',"matrix_H.txt","matrix_L.txt"))
matrixg<-as.data.frame(matri)
colnames(matrixg)<-c('HL','HC1','HC2','dc','LC1','LC2')
matrixg$id<-rownames(matrixg)
forda<-matrixg
write.table(matrixg, "angle_distance.txt", sep="\t")


matrixg<-melt(matrixg,id.vars='id')
pdf('Trajectory_angle_distance.pdf',width=16,height=6)
p<-ggplot(matrixg,aes(x=as.numeric(id),y=value))+geom_line(aes(color=variable,group=variable))+facet_wrap(~variable,ncol=3,scales='free')+xlab('Trajectory')+ylab('Degree or Angstrom')+theme_classic()
print(p)
dev.off()
k<-nrow(dcd);
pdf('Trajectory_angle_distance_histogram.pdf',width=16,height=6)
p<-ggplot(matrixg,aes(x=value,k=k,y=..count..*100/k))+geom_histogram(stat='bin',position='dodge',aes(fill=variable,group=variable),color='gray')+facet_wrap(~variable,ncol=3,scales='free')+ylab('Frequency (%)')+xlab('Degree or Angstrom')+theme_classic()
print(p)
dev.off()

###################Elbow angle########

		matri<-matrix(nrow=nrow(dcd),ncol=1)
		indsH<-atom.select(pdb,'noh',chain='H')
		indsL<-atom.select(pdb,'noh',chain='L')
		hmin=min(pdb$atom[pdb$atom$chain =='H',]$resno)
		lmin=min(pdb$atom[pdb$atom$chain =='L',]$resno)
		matrisc<-matrix(nrow=nrow(dcd),ncol=1)
		inds<-combine.select(indsH,indsL,operator="OR")
		for(i in 1:nrow(dcd)){#
		  write.pdb(pdb=trim(pdb, ind=inds),file='HL.pdb',xyz=dcd[i,inds$xyz],print.segid=T)
		  matri[i]<-system(paste(Traj_folder,'/elbow_angle.pl -f HL.pdb -hend ',seleposeH,' -lend ',seleposeL,sep=''),intern=T)#-hend 121 -lend 323 

		  cat(paste(i,matri[i],'\n',sep="	"))
		}
		system("rm HL.pdb")
		write(matri,file='eblow_angle.txt',sep='')
		da<-data.frame(k=matri)
		pdf('./Trajectory_elbowangle.pdf',width=4,height=4)
		p<-ggplot(da,aes(x=as.numeric(rownames(da)),y=as.numeric(as.character(k))))+geom_line()+xlab('Trajectory')+ylab('Degree')
		print(p)
		dev.off()
		t=nrow(dcd)
		pdf('./Trajectory_eblowangle_histogram.pdf',width=4,height=4)
		p<-ggplot(da,aes(x=as.numeric(as.character(k)),t=t,y=..count..*100/t))+geom_histogram(stat='bin',position='dodge', color='gray',fill='blue')+ylab('Frequency (%)')+xlab('Degree')
		print(p)
		dev.off()

		da$HL<-forda$HL
		pdf('HL_elbow_correlation.pdf',width=8,height=8)
		ggplot(da,aes(x=HL,y=as.numeric(as.character(k))))+geom_point(color='blue',alpha=0.2)+stat_density_2d(aes(fill=..level..),contour=T,geom='polygon')+geom_density2d(color='gray',size=1,alpha=0.2)+scale_color_gradient(low='blue',high='red')+scale_fill_gradient(low='blue',high='red')+theme(axis.text = element_text(size=10))+xlim(-80,-40)+ylim(100,250)+ylab('Elbow angle')
		dev.off()

####interface area and Hbond dynamics#####
#VH-VL bASA
print("VH-VL bASA")
HLa <- atom.select(pdb, chain=c('H','L'))
Hvregion<-atom.select(pdb,chain='H',resno = c(seleposeHst:seleposeH))
Lvregion<-atom.select(pdb,chain='L',resno = c(seleposeLst:seleposeL))
HLa<-combine.select(Hvregion,Lvregion,operator="OR")
pdbHL<-trim(pdb,ind=HLa)
pdbname<-'Traj'
unlink("Fv_Interface_ASA.txt");
for(i in 1:nrow(dcd)){
  fname<-paste('tmp',i,sep='')
  write.pdb(pdb=pdbHL,file=paste(pdbname,'.pdb',sep=''),xyz=dcd[i,HLa$xyz],print.segid=T)
  system(paste(Traj_folder,'/interface_residues.pl -f ',paste(pdbname,'.pdb',sep=''),' -c1 H -c2 L -ASA Fv_Interface_ASA.txt -frame ',i,sep=''))
}
#elbow bASA
Fab <- atom.select(pdb,'noh', chain=c('H','L'))
Fv<-combine.select(Hvregion,Lvregion,operator="OR")
Fvstr<-trim(pdb,ind=Fv)
Fabstr<-trim(pdb,ind=Fab)
Fabstr$atom[Fabstr$atom$resno %in% Fvstr$atom$resno,]$chain<-'H'
Fabstr$atom[!Fabstr$atom$resno %in% Fvstr$atom$resno,]$chain<-'C'

#xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=Fv$xyz,  mobile.inds=ca.inds$xyz)
unlink("Elbow_ASA.txt");
for(i in 1:nrow(dcd)){
  fname<-paste('tmp',i,sep='')
  write.pdb(pdb=Fabstr,file=paste(pdbname,'.pdb',sep=''),xyz=xyz[i,Fab$xyz],print.segid=T)
  system(paste(Traj_folder,'/interface_residues.pl -f ',paste(pdbname,'.pdb',sep=''),' -c1 H -c2 C -ASA Elbow_ASA.txt -frame ',i,sep=''))
}
ASA<-read.table('Elbow_ASA.txt',sep='\t')
colnames(ASA)<-c('Trajectory','chain','Fv_num_res','C_num_res','Fv_SASA','C_SASA','Hbond','Saltbridge')
ASA<-ASA[ASA$Fv_SASA>100,]
pdf('Elbow_ASA_trajectory.pdf',width=8,height=4)
ggplot(ASA[ASA$Fv_SASA>0 & ASA$C_SASA>0,],aes(x=Trajectory))+geom_line(aes(y=Fv_SASA+C_SASA))+ylab('Elbow total SASA')+theme_classic()+theme(axis.text = element_text(size=10))
dev.off()


system(paste(Traj_folder,'/interface_residues_trajectory.pl -f Fv_Interface_ASA.txt -o Interface_Hbond_trajectory.txt', sep=''))
bonds<-read.table('Interface_Hbond_trajectory.txt',head=T)
pdf('Interface_Hbond.pdf',width=12,height=8)
ggplot(bonds,aes(x=trajectory,y=residue_pair))+scale_shape_identity()+geom_point(aes(group=1,color=as.factor(found),shape=124))+scale_color_manual(values=c('white','black'))+facet_wrap(~Bondtype,ncol=1, scales = 'free')+theme_classic()+theme(axis.text = element_text(size=8))
dev.off()

ASA<-read.table('Fv_Interface_ASA.txt',sep='\t')
colnames(ASA)<-c('Trajectory','chain','H_num_res','L_num_res','Heavy_SASA',"Light_SASA",'Hbond','Saltbridge')
ASA<-ASA[ASA$Heavy_SASA>100,]
pdf('Interface_ASA_trajectory.pdf',width=8,height=4)
ggplot(ASA[ASA$Heavy_SASA>0 & ASA$Light_SASA>0,],aes(x=Trajectory))+geom_point(aes(y=Heavy_SASA,color='blue'))+geom_point(aes(y=Light_SASA,color='red'))+scale_color_manual(values=c('blue','red'),label=c('Heavy','Light'))+ylab('SASA')+theme(axis.text = element_text(size=10))
dev.off()

pdf('Interface_ASA_density.pdf',width=8,height=8)
ggplot(ASA[ASA$Heavy_SASA>0 & ASA$Light_SASA>0,],aes(x=Heavy_SASA,y=Light_SASA))+stat_density_2d(aes(fill=..level..),contour=T,geom='polygon')+geom_density2d(aes(),color='gray',size=1,alpha=0.5)+scale_color_gradient(low='blue',high='red')+scale_fill_gradient(low='blue',high='red')+theme(axis.text = element_text(size=10))
dev.off()

allHbond<-read.table('../Hbond_per_pos_total.txt')
allHbond$V2<-factor(allHbond$V2, levels=rev(unique(allHbond$V2)))
pdf('Hbond_per_position.pdf',width=12,height=16)
ggplot(allHbond,aes(x=as.numeric(V1),y=V2))+scale_shape_identity()+geom_point(aes(group=1,color=as.numeric(V3),shape=124))+xlab('Trajectory')+ylab('Position')+scale_color_gradient2(low='white',high='red')+theme(axis.text.y = element_text(size=3))
dev.off()

unlink(c('change.pdb','change_interface_residue.txt','Traj.pdb'))

quit()

