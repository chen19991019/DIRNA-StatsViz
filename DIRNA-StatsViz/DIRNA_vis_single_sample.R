library(optparse)
library(ggplot2)
library(hexbin)
library(ggpubr)
library(ggseqlogo)
make_heatmap = function(sample,segment)
{
  bw = 50
  file = paste(sample,"/",sample,"_",segment,"_dels_map_density.txt",sep="")
  data = read.csv(file, sep="\t", header=TRUE)
  ggp = ggplot(data, aes(x=DelStart, y=DelEnd) ) +
    geom_bin2d(binwidth=bw) +
    scale_fill_gradient(low = "yellow", high = "red", trans="log10") +
    theme_minimal() +
    labs(title = segment)+
    xlab("Deletion Start") +
    ylab("Deletion End")
  return(ggp)
}
make_bar_for_breakpoints = function(Nucleotide)
{
  master_e=as.data.frame(t(data_e[Nucleotide,2:5])) #等是5个的时候改成2:11
  master_e$Type="Expected"
  master_e$Position=row.names(master_e)
  master_o=as.data.frame(t(data_o[Nucleotide,2:5])) #等是5个的时候改成2:11
  master_o$Type="Observed"
  master_o$Position=row.names(master_o)
  master=rbind(master_e,master_o)
  colnames(master)[1]="Percentage"
  master$Position = factor(master$Position,levels=c("up_2","up_1","down_1","down_2"))
  #5个的时候"up_5","up_4","up_3","up_2","up_1","down_1","down_2","down_3","down_4","down_5"
  ggp=ggplot(master,aes(Position,Percentage,fill=Type))+
    geom_bar(stat="identity",position="dodge")+
    scale_fill_manual(values = c('#bdbdbd','#fe9929'))+
    theme_bw()+
    labs(title=paste("Nucleotide_",Nucleotide,sep=""))+
    ylab("Nucleotide Occurrence(%)")+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))
  return(ggp)
}

parameter=commandArgs(trailingOnly = TRUE)
sample=parameter[1]
workdirectory=parameter[2]



#When performing analysis between projects, the path to each sample folder will change, 
#so the main path will need to be changed.
mainpath=getwd()
if (is.na(workdirectory)){
  workdirectory=mainpath
}else{
  workdirectory=paste(mainpath,"/",workdirectory,sep="")
}
setwd(workdirectory)

####dirna length plot####
pb2=read.csv(paste(sample,"/dirna_length/PB2_dirna_length.csv",sep=""),header=TRUE,sep= ",")
pb1=read.csv(paste(sample,"/dirna_length/PB1_dirna_length.csv",sep=""),header=TRUE,sep= ",")
pa=read.csv(paste(sample,"/dirna_length/PA_dirna_length.csv",sep=""),header=TRUE,sep= ",")
ha=read.csv(paste(sample,"/dirna_length/HA_dirna_length.csv",sep=""),header=TRUE,sep= ",")
np=read.csv(paste(sample,"/dirna_length/NP_dirna_length.csv",sep=""),header=TRUE,sep= ",")
na=read.csv(paste(sample,"/dirna_length/NA_dirna_length.csv",sep=""),header=TRUE,sep= ",")
m=read.csv(paste(sample,"/dirna_length/M_dirna_length.csv",sep=""),header=TRUE,sep= ",")
ns=read.csv(paste(sample,"/dirna_length/NS_dirna_length.csv",sep=""),header=TRUE,sep= ",")

acname=paste(sample,"_coverage.txt", sep="")
ac=read.csv(paste(sample,"/",acname, sep=""),header=TRUE,sep= ",")
ac[6,1]="NA"
rownames(ac)=ac$RefName
ha$Name="HA"
pa$Name="PA"
pb1$Name="PB1"
pb2$Name="PB2"
m$Name="M"
na$Name="NA"
np$Name="NP"
ns$Name="NS"
master=rbind(ha,pa,pb1,pb2,m,na,np,ns)
master$Percentage=master$Percentage*100
master$Name=factor(master$Name,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
master$Length=factor(master$Length,levels=c("1<=length<50","50<=length<100","100<=length<250","250<=length<500","500<=length<1000","1000<=length<1500","1500<=length<2000","2000<=length"))
length_plot=ggplot(master, aes(x = Name, y = Percentage, fill = Length)) +  # Create stacked bar chart (删掉了label = Percentage)
  geom_bar(stat = "identity", width = 0.7)+
  scale_fill_brewer(palette = "YlGnBu")+ 
  theme_bw()+
  labs(x='Segment',title = "DIRNA_length_frequency_plot")
length_plot
ggsave(length_plot,filename = paste(sample,"/plots/length_frequency_plot.pdf",sep=""),height = 6, width = 8)

####Count and Scaled####
count_scale=data.frame(Name=c(rep("PB2",2),rep("PB1",2),rep("PA",2),rep("HA",2),rep("NP",2),rep("NA",2),rep("M",2),rep("NS",2)),
                       Type=c("Count","Scaled","Count","Scaled","Count","Scaled","Count","Scaled","Count","Scaled","Count","Scaled","Count","Scaled","Count","Scaled"),
                       Number=c(sum(pb2$Count),round((sum(pb2$Count)/ac["PB2","AvDepth"])*5000),
                                sum(pb1$Count),round((sum(pb1$Count)/ac["PB1","AvDepth"])*5000),
                                sum(pa$Count),round((sum(pa$Count)/ac["PA","AvDepth"])*5000),
                                sum(ha$Count),round((sum(ha$Count)/ac["HA","AvDepth"])*5000),
                                sum(np$Count),round((sum(np$Count)/ac["NP","AvDepth"])*5000),
                                sum(na$Count),round((sum(na$Count)/ac["NA","AvDepth"])*5000),
                                sum(m$Count),round((sum(m$Count)/ac["M","AvDepth"])*5000),
                                sum(ns$Count),round((sum(ns$Count)/ac["NS","AvDepth"])*5000)))
count_scale$Name=factor(count_scale$Name,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
count_scale_plot=ggplot(count_scale,aes(Name,Number,fill=Type))+
  geom_bar(stat="identity",position = 'dodge',width = 0.8)+
  scale_fill_manual(values = c('#c6dbef','#4292c6'))+
  theme_bw()+
  scale_y_continuous(name = "Count",
                     sec.axis = sec_axis(~./5000,name = "Scaled"))+
  labs(x='Segment',y='Frequency',title="DIRNA_Count_Scale")
count_scale_plot
ggsave(count_scale_plot,filename = paste(sample,"/plots/count_scale_plot.pdf",sep=""),height = 6, width = 8)

####dirna length percentage plot####

pb2=read.csv(paste(sample,"/dirna_length_percentage/PB2_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
pb1=read.csv(paste(sample,"/dirna_length_percentage/PB1_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
pa=read.csv(paste(sample,"/dirna_length_percentage/PA_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
ha=read.csv(paste(sample,"/dirna_length_percentage/HA_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
np=read.csv(paste(sample,"/dirna_length_percentage/NP_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
na=read.csv(paste(sample,"/dirna_length_percentage/NA_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
m=read.csv(paste(sample,"/dirna_length_percentage/M_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")
ns=read.csv(paste(sample,"/dirna_length_percentage/NS_dirna_length_percentage.csv",sep=""),header=TRUE,sep= ",")


ha$Name="HA"
pa$Name="PA"
pb1$Name="PB1"
pb2$Name="PB2"
m$Name="M"
na$Name="NA"
np$Name="NP"
ns$Name="NS"
master=rbind(ha,pa,pb1,pb2,m,na,np,ns)
master$Percentage=master$Percentage*100
master$Name=factor(master$Name,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
master$Length_Percentage=factor(master$Length_Percentage,levels=c("0<=length%<20","20<=length%<40","40<=length%<50","50<=length%<60","60<=length%<70","70<=length%<80","80<=length%<90","90<=length%<=100"))
length_percentage_plot=ggplot(master, aes(x = Name, y = Percentage, fill = Length_Percentage)) +  # Create stacked bar chart
  geom_bar(stat = "identity", width = 0.7)+
  scale_fill_brewer(palette = "YlGnBu")+
  theme_bw()+
  labs(x='Segment',y='Frequency',title = "DIRNA_length_percentage_frequency_plot")
length_percentage_plot
ggsave(length_percentage_plot,filename = paste(sample,"/plots/length_percentage_frequency_plot.pdf",sep=""),height = 6, width = 8)

####for remaining length####

pb2=read.csv(paste(sample,"/remaining_length/PB2_remaining_length.csv",sep=""),header=TRUE,sep= ",")
pb1=read.csv(paste(sample,"/remaining_length/PB1_remaining_length.csv",sep=""),header=TRUE,sep= ",")
pa=read.csv(paste(sample,"/remaining_length/PA_remaining_length.csv",sep=""),header=TRUE,sep= ",")
ha=read.csv(paste(sample,"/remaining_length/HA_remaining_length.csv",sep=""),header=TRUE,sep= ",")
np=read.csv(paste(sample,"/remaining_length/NP_remaining_length.csv",sep=""),header=TRUE,sep= ",")
na=read.csv(paste(sample,"/remaining_length/NA_remaining_length.csv",sep=""),header=TRUE,sep= ",")
m=read.csv(paste(sample,"/remaining_length/M_remaining_length.csv",sep=""),header=TRUE,sep= ",")
ns=read.csv(paste(sample,"/remaining_length/NS_remaining_length.csv",sep=""),header=TRUE,sep= ",")

ha$Name="HA"
pa$Name="PA"
pb1$Name="PB1"
pb2$Name="PB2"
m$Name="M"
na$Name="NA"
np$Name="NP"
ns$Name="NS"
master=rbind(ha,pa,pb1,pb2,m,na,np,ns)
master$Percentage=master$Percentage*100
master$Name=factor(master$Name,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
master$Remaining_Length=factor(master$Remaining_Length,levels=c("1<=length<50","50<=length<100","100<=length<250","250<=length<500","500<=length<1000","1000<=length<1500","1500<=length<2000","2000<=length"))
remaining_length_plot=ggplot(master, aes(x = Name, y = Percentage, fill = Remaining_Length)) +  # Create stacked bar chart
  geom_bar(stat = "identity", width = 0.7)+
  scale_fill_brewer(palette = "YlGnBu")+
  theme_bw()+
  labs(x='Segment',y='Frequency',title = "Remaining_length_frequency_plot")
remaining_length_plot
ggsave(remaining_length_plot,filename = paste(sample,"/plots/remaining_length_frequency_plot.pdf",sep=""),height = 6, width = 8)

####combining_percentage bar plot####
combining_percentage_barplot=ggarrange(length_plot, length_percentage_plot, remaining_length_plot, 
                                      ncol = 3, nrow = 1)
ggsave(combining_percentage_barplot,filename = paste(sample,"/plots/combining_frequency_barplot.pdf", sep=""), height = 8, width = 20)

####heatmap####
for(se in c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
{
  ggp = make_heatmap(sample,se)
  ggsave(ggp,filename = paste(sample,"/plots/heatmap_",se,".pdf", sep=""), height = 6, width = 7)
}

####combining_heatmap####
pb2plot = make_heatmap(sample,"PB2")
pb1plot = make_heatmap(sample,"PB1")
paplot = make_heatmap(sample,"PA")
haplot = make_heatmap(sample,"HA")
naplot = make_heatmap(sample,"NA")
nsplot = make_heatmap(sample,"NS")
mplot = make_heatmap(sample,"M")
npplot = make_heatmap(sample,"NP")
combining_heatmap=ggarrange(pb2plot, pb1plot, paplot, haplot, npplot, naplot, mplot, nsplot, 
                  ncol = 2, nrow = 4)
ggsave(combining_heatmap,filename = paste(sample,"/plots/combining_heatmap.pdf", sep=""), height = 12, width = 8)
####breakpoints####
data_e=read.csv(paste(sample,"/breakpoints/breakpoints_expected_percentage.csv",sep=""),header=TRUE,sep=",")
row.names(data_e)=data_e$Nucleotide
data_o=read.csv(paste(sample,"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
row.names(data_o)=data_e$Nucleotide

a=make_bar_for_breakpoints("A")
g=make_bar_for_breakpoints("G")

c=make_bar_for_breakpoints("C")
t=make_bar_for_breakpoints("T")

combining_breakpoint=ggarrange(a,g,c,t, ncol = 2, nrow = 2)
ggsave(combining_breakpoint,filename = paste(sample,"/plots/combining_breakpoints.pdf",sep=""),height = 6, width = 8)

####seq WebLogo(frequency of each amino or nucleic acid)####
seqdata=read.table(paste(sample,"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=FALSE,sep=",")
seqname=seqdata[2:5,1]
seqdata <- rbind(
  as.numeric(seqdata[2,2:5]),
  as.numeric(seqdata[3,2:5]),
  as.numeric(seqdata[4,2:5]),
  as.numeric(seqdata[5,2:5])
)
row.names(seqdata)=seqname
seq=ggseqlogo(seqdata, method="bits")+
  xlab("Position")+
  theme(axis.text.x = element_blank())+
  labs(title="Frequency_of_nucleotides(A/G/C/T)")+
  annotate("text", x = 1, y = -0.03, label = "up_2",color="black",size = 4)+
  annotate("text", x = 2, y = -0.03, label = "up_1",color="black",size = 4)+
  annotate("text", x = 3, y = -0.03, label = "down_1",color="black",size = 4)+
  annotate("text", x = 4, y = -0.03, label = "down_2",color="black",size = 4)+
  guides(colour = "none")
seq
ggsave(seq,filename = paste(sample,"/plots/frequency_of_nucleotides.pdf",sep=""),height = 6, width = 8)


#R --vanilla --slave --args minion < /Users/chenhong1/DIRNA-StatsViz/DIRNA_vis_single_sample.R


