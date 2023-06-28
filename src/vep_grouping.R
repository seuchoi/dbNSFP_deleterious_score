vep_grouping<-function(infile,types=c("hclof","hclof_noflag","missense","synonymous"),dbnsfp=c(4.2,4.3,4.4),mintools=7){

## infile = VEP annotated file using LOFTEE and dbnsfp
## types = high confidence loss of function variants (hclof), high confidence loss of function variants without flag (hclof_noflag), missense variant, synonymous variants
## dbnsfp = we are supporting dbnsfp 4.2, 4.3, 4.4 versions
## mintools = the minimun number of bioinformatics tools to incldue in the grouping file


header0<-fread(cmd=paste0("zcat ",infile," | head -n 500 | grep ^#Uploaded_variation "),header=F,data.table=F,sep="\t")

basefilename<-basename(infile)

if(any(c("hclof","hclof_noflag") %in% types)){

dat0<-fread(cmd=paste0("zcat ",infile," | egrep HC "),header=F,data.table=F,sep="\t")

names(dat0)<-header0[1,]

#####
##### select protein coding genes
colns<-names(dat0)
avail<-length(which(colns=="BIOTYPE"))

if(!isZero(avail)){
dat0<-subset(dat0,BIOTYPE=="protein_coding")
}else{
message(" BIOTYPE variable is not available. Unable to select protein coding genes.. ")
}


#####
#####
if("hclof" %in% types){
sub1<-subset(dat0,LoF=="HC")
score0<-try(delscore_af_gnomad_data_dbnsfp(sub1,version=dbnsfp,delvar="Dprop",mintool=mintools,in.transcript=TRUE),silent=T)
group<-score0
names(group)[5]<-"group_id"
##### save file
outfile<-paste0(basefilename,".hclof.RData")
save(group,file=outfile)

}

if("hclof_noflag" %in% types){
sub1<-subset(dat0,LoF=="HC" & LoF_flags=="-" )
score0<-try(delscore_af_gnomad_data_dbnsfp(sub1,version=dbnsfp,delvar="Dprop",mintool=mintools,in.transcript=TRUE),silent=T)
group<-score0
names(group)[5]<-"group_id"
##### save file
outfile<-paste0(basefilename,".hclof_noflag.RData")
save(group,file=outfile)
}
}

if("missense" %in% types){

dat0<-fread(cmd=paste0("zcat ",infile," | egrep missense "),header=F,data.table=F,sep="\t")
names(dat0)<-header0[1,]

#####
##### select protein coding genes
colns<-names(dat0)
avail<-length(which(colns=="BIOTYPE"))

if(!isZero(avail)){
dat0<-subset(dat0,BIOTYPE=="protein_coding")
}else{
message(" BIOTYPE variable is not available. Unable to select protein coding genes.. ")
}

#####
#####
sub1<-dat0[grep("missense_variant",dat0$Consequence),]
score0<-try(delscore_af_gnomad_data_dbnsfp(sub1,version=dbnsfp,delvar="Dprop",mintool=mintools,in.transcript=TRUE),silent=T)
group<-score0
names(group)[5]<-"group_id"
##### save file
outfile<-paste0(basefilename,".missense.RData")
save(group,file=outfile)
}


if("synonymous" %in% types){

dat0<-fread(cmd=paste0("zcat ",infile," | egrep synonymous "),header=F,data.table=F,sep="\t")
names(dat0)<-header0[1,]

#####
##### select protein coding genes
colns<-names(dat0)
avail<-length(which(colns=="BIOTYPE"))

if(!isZero(avail)){
dat0<-subset(dat0,BIOTYPE=="protein_coding")
}else{
message(" BIOTYPE variable is not available. Unable to select protein coding genes.. ")
}

#####
#####
sub1<-dat0[grep("synonymous_variant",dat0$Consequence),]
score0<-try(delscore_af_gnomad_data_dbnsfp(sub1,version=dbnsfp,delvar="Dprop",mintool=mintools,in.transcript=TRUE),silent=T)
group<-score0
names(group)[5]<-"group_id"
##### save file
outfile<-paste0(basefilename,".synonymous.RData")
save(group,file=outfile)
}

}
