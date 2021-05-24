SSNChen=function(ds,fn="SSNChen.png",width = 18, height =5,dpi = 600,
descr="FASN (Moschus berezovskii)",Author="Drafted by Chen 2021"){
address=getwd() # Previous record location

library(ggplot2) # Loading ggplot2

cds=read.csv(deparse(substitute(ds))) # Reading sheet
# Data structure included that R_Lengths(number),R_Names(string),starCDS(number),endCDS(number),PloyASignal(number),KeySites(number),miRTarNames(string),miRTarSites(number),Psites(number),NPrimer(string),FPrimer(number),RPrimer(number)

omit=function(cc){a=cc[!is.na(cc)];a[which(a!="")]} # NA omit in sheet

AP=function(org,NP,FP,RP,c=1.28,d=1.32){
n=rep(c(0,0.1),ceiling(length(NP)/2))[1:length(NP)]
org+geom_rect(aes(xmin=FP,xmax=RP,ymin=c+n,ymax=d+n),fill="grey")+
geom_text(aes(x=(FP+RP)/2,y=c+n+0.025,label=NP))} # Primers adding

AM=function(org,mirs,labs,a=0.92,b=-0.05,s=2){ 
org+geom_text(aes(x=mirs,y=rep(a,length(mirs))+seq(0,b*(length(mirs)-1),b),label=labs,vjust="left"),fontface="bold",color="blue1",
size=s,inherit.aes =F,check_overlap = T)} # Adding miRNA information

# Loading data from sheet
x=omit(cds$R_Lengths)        # length of each setting regions
y="Nucleic acid"             # fixed sequence name NO USE
zt=omit(cds$R_Names)         # names of each regions
starCDS=omit(cds$starCDS)    # CDS star site
endCDS=omit(cds$endCDS)      # CDS end site
AATAAA=omit(cds$PloyASignal) # AATAAA PloyA signal
mp=omit(cds$KeySites)        # key sites you need
mirsites=omit(cds$miRTarSites)  # miRNA targetted sites
mirlabs=paste("¡ô",omit(cds$miRTarNames),sep="")  #miRNA names
PS=omit(cds$Psites)   # Phosphorylation Sites you need
NP=omit(cds$NPrimer)  # Primer names
FP=omit(cds$FPrimer)  # Forward Primer first site
RP=omit(cds$RPrimer)  # Reverse Primer last site
IP=omit(cds$InterP)   # Interaction proteins

# Ready information
linecol="black"    # Border color
zs=paste("S",paste(seq(1,length(x)),zt,sep=": "),sep="") # Adding signal names
z=paste(toupper(letters)[rev(seq(1,length(x)))],zt,sep=":") # Sort annotation by letters
inf=data.frame(Length=x,NASeq=y,Section=z) # Simple frame generation
bcolor=c("#CDBE70","#CD6839","#CD2626","#FFB5C5","#7D26CD",
"#66CDAA","#CD8500","#CDCD00","#00CD00","#668B8B",
"#9ACD32","#8B658B","#CD00CD")[1:length(x)] # Candidate colors

# Draft core orgin
orgin=with(inf,ggplot()+theme(axis.text.y=element_blank(),
panel.background = element_blank(),panel.border = element_blank(),
text=element_text(face="bold"),axis.text = element_text(face="bold"),
legend.text=element_text(face="bold"),legend.key=element_rect(colour='black',size=1),
legend.position="top")+
scale_fill_manual(values =bcolor,labels=rev(zs))+
guides(fill=guide_legend(reverse=TRUE))+
geom_bar(aes(Length,NASeq,fill=Section),stat="identity",position = "stack",width=0.1)+
geom_point(aes(x=mp,y=1),size=3, shape=10, color="black",show.legend=F, inherit.aes=F)+
geom_rect(aes(xmin=starCDS,xmax=endCDS,ymin=1.08,ymax=1.12),fill="DeepSkyBlue2")+
geom_text(aes(x= (endCDS-starCDS)/2+starCDS,y=1.12,
label = "CDS Region",vjust=-0.25),fontface="bold")+
geom_text(aes(x= AATAAA,y=1.05,label = "¡ýPloyA",vjust=-0.5))+
geom_line(aes(x=c(0,sum(Length)),y=0.95),size=1,color=linecol)+
geom_line(aes(x=c(0,sum(Length)),y=1.05),size=1,color=linecol)+
geom_line(aes(x=0,y=c(0.95,1.05)),size=1,color=linecol)+
geom_line(aes(x=sum(Length),y=c(0.95,1.05)),size=1,color=linecol)+
geom_text(aes(x=sum(Length)/2,y=0.75,label=descr),fontface="bold.italic",size=5)+
geom_point(aes(x=PS,y=1.105),size=3, shape=18, color="orange",
show.legend=F, inherit.aes =F)+
geom_point(aes(x=5,y=0.85),size=3, shape=10, color="black",inherit.aes =F)+
geom_text(aes(x=5,y=0.85,label="Key Sites",hjust=-0.13),fontface="bold")+
geom_point(aes(x=5,y=0.75),size=3, shape=18, color="orange",inherit.aes =F)+
geom_text(aes(x=5,y=0.75,label="Phosphorylation Sites",hjust=-0.05),fontface="bold")+
geom_text(aes(x=5,y=0.62,label=paste("InterPs: ",paste(IP,collapse="; ")),hjust=-0.001),fontface="bold",size=3)+
geom_text(aes(x=sum(Length)*0.95,y=0.45,label=Author)))

p=AP(orgin,NP,FP,RP) # Adding primer mode
p=AM(p,mirsites,mirlabs) # Adding miRNA information

# Save files
if(file.exists("./SSNDraft")==TRUE){cat("SSNDraft Existed!\n")}else{
dir.create("./SSNDraft",recursive=TRUE);cat("SSNDraft Created\n")}
setwd("./SSNDraft")
ggsave(paste(gsub(":","_",Sys.time()),fn,sep="_"),p,width = width,height=height,dpi=dpi)
cat("File had been saved sucessfully under ",paste(address,"/SSNDraft",sep=""),"\n")
setwd(address);p}