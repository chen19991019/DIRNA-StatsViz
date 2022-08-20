library(optparse)
library(ggplot2)
library(dplyr)
library(grid)
library(UpSetR)
library(ggvenn)
library(ggpubr)
library(ggrepel)
parameter=commandArgs(trailingOnly = TRUE)
sampleconfigure=parameter[1]
upset_threshold=venn_threshold=parameter[2]
workdirectory=parameter[3]
####Function####
make_venn_boxplot = function(segment)
{
  y=list()
  for(i in co$Samplename)
  {
    filepath=paste(i,"/",i,"_",segment,"_dels_map.txt",sep="")
    i_data=read.table(filepath,header=TRUE,sep="\t")
    x=list(
      i_list=i_data[which(i_data$Count>venn_threshold),]$DelStart.DelEnd
    )
    names(x)=i
    y=append(y,x)
  }
  
  ggp=ggvenn(y, stroke_size = 0.5, set_name_size = 4)+
    labs(title = segment)
  return(ggp)
}

make_upset_boxplot = function(segment)
{
  t=number_of_sample*0.5
  y=list()
  for(i in co$Samplename)
  {
    filepath=paste(i,"/",i,"_",segment,"_dels_map.txt",sep="")
    i_data=read.table(filepath,header=TRUE,sep="\t")
    x=list(
      i_list=i_data[which(i_data$Count>upset_threshold),]$DelStart.DelEnd
    )
    if (length(unlist(x)) == 0){
      next
    }else{
      names(x)=i
      y=append(y,x)
    }
  }
  if (length(y) >= 2){
    y=fromList(y)
    ggp=upset(y,
              nsets=length(y),
              nintersects = 1000,
              keep.order = TRUE,
              sets.x.label = "Total Type",
              mainbar.y.label = paste("Type of DiRNA for ",segment,sep=""),
              point.size = 6.5-t,
              line.size = 4.5-t,
              number.angles = 0,
              text.scale = c(5-t, 4.7-t, 4.7-t, 4.5-t, 5-t, 4.5-t),
              # = c(3.5, 3.2, 3.2, 3, 3.5, 3),
              order.by="degree",
              matrix.color='#4285F4',
              main.bar.color = 'black',
              sets.bar.color = "#bdbdbd",
              set_size.show=TRUE)
    return(ggp)
  }else{
    ggp = paste("!!!Warning:There is one sample where the count of all DiRNA of the ",segment," segment is less than the threshold, so the upset plot does not apply to segment ",segment,sep="")
    return(ggp)
  }
}

make_line_plot_for_Count = function(number_of_sample,segment)
{
  summary=read.table("Summary.csv",header=TRUE,sep=",")
  row.names(summary)=summary$Segment
  sumend=9+number_of_sample
  sum_count=summary[segment,10:sumend]
  sum_count=as.data.frame(t(sum_count))
  sum_count$Samplename=co$Samplename
  sum_count$Type="All_DiRNA"
  
  sharedstart=sumend+1
  sharedend=9+number_of_sample+number_of_sample
  shared_count=summary[segment,sharedstart:sharedend]
  shared_count=as.data.frame(t(shared_count))
  shared_count$Samplename=co$Samplename
  shared_count$Type="Shared_DiRNA"
  master=rbind(sum_count,shared_count)
  colnames(master)[1]="Counts"
  master$Samplename=factor(master$Samplename,levels=c(co$Samplename))
  
  ggp=ggplot(master,aes(Samplename,Counts,color=Type,group=Type,label = round(Counts, 1)))+
    geom_point(size=2)+
    geom_line(cex=1)+
    geom_text_repel()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
    labs(title = paste(segment,"_Count",sep=""))+
    scale_color_manual(values = c('#a50f15','#08519c'))
  return(ggp)
}

make_line_plot_for_Type = function(number_of_sample,segment)
{
  summary=read.table("Summary.csv",header=TRUE,sep=",")
  row.names(summary)=summary$Segment
  sum_star=9+number_of_sample+number_of_sample+1
  sumend=9+number_of_sample+number_of_sample+number_of_sample
  sum_dirna=summary[segment,sum_star:sumend]
  sum_dirna=as.data.frame(t(sum_dirna))
  sum_dirna$Samplename=co$Samplename
  sum_dirna$Type="All_DiRNA"
  colnames(sum_dirna)[1]="Types"
  
  shared_dirna=as.data.frame(matrix(summary[segment,5], ncol = 1, nrow = nrow(sum_dirna)))
  shared_dirna$Samplename=co$Samplename
  shared_dirna$Type="Shared_DiRNA"
  colnames(shared_dirna)[1]="Types"
  master=rbind(sum_dirna,shared_dirna)
  master$Samplename=factor(master$Samplename,levels=c(co$Samplename))
  
  ggp=ggplot(master,aes(Samplename,Types,color=Type,group=Type,label = round(Types, 1)))+
    geom_point(size=2)+
    geom_line(cex=1)+
    geom_text_repel()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1))+
    labs(title = paste(segment,"_DiRNA",sep=""))+
    scale_color_manual(values = c('#a50f15','#08519c'))
  return(ggp)
}


make_dot_plot_for_breakpoints = function(Nucleotide)
{
  expected=read.table("breakpoints/breakpoints_expected_percentage.csv",header=TRUE,sep=",")
  row.names(expected)=expected$Nucleotide
  expected_data=as.data.frame(t(expected[Nucleotide,2:5]))
  expected_data$Type="Expected"
  expected_data$Position=row.names(expected_data)
  
  sample1=read.table(paste(co[1,1],"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
  row.names(sample1)=sample1$Nucleotide
  sample1_data=as.data.frame(t(sample1[Nucleotide,2:5]))
  sample1_data$Type=co[1,1]
  sample1_data$Position=row.names(sample1_data)
  master=rbind(expected_data,sample1_data)
  for(samplename in co[2:number_of_sample,1])
  {
    sample=read.table(paste(samplename,"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
    row.names(sample)=sample$Nucleotide
    sample_data=as.data.frame(t(sample[Nucleotide,2:5]))
    sample_data$Type=samplename
    sample_data$Position=row.names(sample_data)
    master=rbind(master,sample_data)
  }
  colnames(master)[1]="Percentage"
  master$Position =factor(master$Position,levels=c("up_2","up_1","down_1","down_2"))
  ggp=ggplot(master, aes(x=Position, y=Percentage, col=Type)) +
    geom_point(size = 4, alpha = 0.5)+
    theme_bw()+
    labs(title=paste("Nucleotide_",Nucleotide,sep=""))+
    ylab("Nucleotide Occurrence(%)")
  return(ggp)
}

save_plot = function(plot,path_to_plot_folder_plotname,suffix) 
{
  pdf(paste(path_to_plot_folder_plotname,suffix,".pdf", sep=""), height = 8, width = 20)
  print(plot)
  dev.off()
}
####set the workdirectory####
mainpath=getwd()
if (is.na(workdirectory)){
  workdirectory=mainpath
}else{
  workdirectory=paste(mainpath,"/",workdirectory,sep="")
}
setwd(workdirectory)
####read configure file and get the number of samples####
co=read.csv(sampleconfigure,header=TRUE,sep= ",")
number_of_sample=length(co$Samplename)
####summary bar####
summary=read.csv("Summary.csv",header=TRUE,sep= ",")
sum=as.data.frame(summary$Sum_Counts)
colnames(sum)[1]="Counts"
sum$Segment=summary$Segment
sum$Type="All_DiRNA"
shared=as.data.frame(summary$Shared_Counts*100)
colnames(shared)[1]="Counts"
shared$Segment=summary$Segment
shared$Type="Shared_DiRNA"
master=rbind(sum,shared)
master$Segment=factor(master$Segment,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
master$Type=factor(master$Type,levels=c("All_DiRNA","Shared_DiRNA"))
summary_count=ggplot(master,aes(Segment,Counts,fill=Type))+
  geom_bar(stat="identity",position = 'dodge',width = 0.8)+
  scale_fill_manual(values = c('#2171b5','#b10026'))+
  theme_bw()+
  scale_y_continuous(name = "All_Counts",
                     sec.axis = sec_axis(~./100,name = "Shared_Counts"))+
  labs(x='Segment',title="Counts_of_Sum_and_Shared_DiRNA")
ggsave(summary_count,filename = "plots/summary_count.pdf",height = 6, width = 8)

summary=read.csv("Summary.csv",header=TRUE,sep= ",")
sum=as.data.frame(summary$Sum_Types)
colnames(sum)[1]="Types"
sum$Segment=summary$Segment
sum$Type="All_DiRNA"
shared=as.data.frame(summary$Shared_Types*100)
colnames(shared)[1]="Types"
shared$Segment=summary$Segment
shared$Type="Shared_DiRNA"

master=rbind(sum,shared)
master$Segment=factor(master$Segment,levels=c("PB2","PB1","PA","HA","NP","NA","M","NS"))
master$Type=factor(master$Type,levels=c("All_DiRNA","Shared_DiRNA"))
summary_type=ggplot(master,aes(Segment,Types,fill=Type))+
  geom_bar(stat="identity",position = 'dodge',width = 0.8)+
  scale_fill_manual(values = c('#2171b5','#b10026'))+
  theme_bw()+
  scale_y_continuous(name = "All_Types",
                     sec.axis = sec_axis(~./100,name = "Shared_Types"))+
  labs(x='Segment',title="Types_of_All_and_Shared_DiRNA")
ggsave(summary_type,filename = "plots/summary_type.pdf",height = 6, width = 8)
combining_summary=ggarrange(summary_count, summary_type,labels=c("A","B"), ncol = 2, nrow = 1)
ggsave(combining_summary,filename = "plots/combining_summary.pdf",height = 6, width = 14)
####venn####
#combining_venn#
if (number_of_sample<5){
  for(se in c("PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"))
  {
    ggp = make_venn_boxplot(se)
    ggsave(ggp,filename = paste("plots/venn_",se,".pdf",sep=""),height = 6, width = 8)
  }
  pb2venn = make_venn_boxplot("PB2")
  pb1venn = make_venn_boxplot("PB1")
  pavenn = make_venn_boxplot("PA")
  havenn = make_venn_boxplot("HA")
  npvenn = make_venn_boxplot("NP")
  navenn = make_venn_boxplot("NA")
  mvenn = make_venn_boxplot("M")
  nsvenn = make_venn_boxplot("NS")
  combining_venn=ggarrange(pb2venn, pb1venn, pavenn, havenn, npvenn, navenn, mvenn, nsvenn, ncol = 2, nrow = 4)
  ggsave(combining_venn,filename = "plots/combining_venn.pdf",height = 16, width = 12)
  
}else{
  print("!!!Warning:The number of sample is more than 4, venn plot cannot be used to visualise the shared DiRNA between samples.")
}


####upset####
if (number_of_sample>7){
  print("!!!Warning:The number of sample is more than 7, upset plot cannot be used to visualise the shared DiRNA between samples.")
}else{
  for(se in c("NP","PB2", "PB1", "PA", "HA", "NA", "M", "NS"))
  {
    ggp = make_upset_boxplot(se)
    if (is.character(ggp)){
      print(ggp)
    }else{
      save_plot(ggp,"plots/upset_",se)
    }
  }
}
#when there is not  shared dirna within samples, 
#it will output that"geom_path: Each group consists of only one observation. Do you need to adjust the group aesthetic?".

####line_for_Count####
pb2 = make_line_plot_for_Count(number_of_sample,"PB2")
ggsave(pb2,filename = "plots/line_count_PB2.pdf",height = 6, width = 8)

pb1 = make_line_plot_for_Count(number_of_sample,"PB1")
ggsave(pb1,filename = "plots/line_count_PB1.pdf",height = 6, width = 8)

pa = make_line_plot_for_Count(number_of_sample,"PA")
ggsave(pa,filename = "plots/line_count_PA.pdf",height = 6, width = 8)

ha = make_line_plot_for_Count(number_of_sample,"HA")
ggsave(ha,filename = "plots/line_count_HA.pdf",height = 6, width = 8)

np = make_line_plot_for_Count(number_of_sample,"NP")
ggsave(np,filename = "plots/line_count_NP.pdf",height = 6, width = 8)

na = make_line_plot_for_Count(number_of_sample,"Na")  
ggsave(na,filename = "plots/line_count_Na.pdf",height = 6, width = 8)

m = make_line_plot_for_Count(number_of_sample,"M")
ggsave(m,filename ="plots/line_count_M.pdf",height = 6, width = 8)

ns = make_line_plot_for_Count(number_of_sample,"NS")
ggsave(ns,filename = "plots/line_count_NS.pdf",height = 6, width = 8)

#combining_line#
combining_line_for_Count=ggarrange(pb2, pb1, pa, ha, np, na, m, ns, ncol = 2, nrow = 4)
ggsave(combining_line_for_Count,filename = "plots/combining_line_for_DiRNA_count.pdf", height = 18, width = 12)

####line_for_Type####
pb2 = make_line_plot_for_Type(number_of_sample,"PB2")
ggsave(pb2,filename = "plots/line_DiRNA_PB2.pdf",height = 6, width = 8)

pb1 = make_line_plot_for_Type(number_of_sample,"PB1")
ggsave(pb1,filename = paste("plots/line_DiRNA_","PB1",".pdf",sep=""),height = 6, width = 8)

pa = make_line_plot_for_Type(number_of_sample,"PA")
ggsave(pa,filename = paste("plots/line_DiRNA_","PA",".pdf",sep=""),height = 6, width = 8)

ha = make_line_plot_for_Type(number_of_sample,"HA")
ggsave(ha,filename = paste("plots/line_DiRNA_","HA",".pdf",sep=""),height = 6, width = 8)

np = make_line_plot_for_Type(number_of_sample,"NP")
ggsave(np,filename = paste("plots/line_DiRNA_","NP",".pdf",sep=""),height = 6, width = 8)

na = make_line_plot_for_Type(number_of_sample,"Na")
ggsave(na,filename = paste("plots/line_DiRNA_","NA",".pdf",sep=""),height = 6, width = 8)

m = make_line_plot_for_Type(number_of_sample,"M")
ggsave(m,filename = paste("plots/line_DiRNA_","M",".pdf",sep=""),height = 6, width = 8)

ns = make_line_plot_for_Type(number_of_sample,"NS")
ggsave(ns,filename = paste("plots/line_DiRNA_","NS",".pdf",sep=""),height = 6, width = 8)
#combining_line#
combining_line_for_Type=ggarrange(pb2, pb1, pa, ha, np, na, m, ns, ncol = 2, nrow = 4)
ggsave(combining_line_for_Type,filename = "plots/combining_line_for_DiRNA_type.pdf", height = 18, width = 12)

####dot plot for multiple smaples####
a=make_dot_plot_for_breakpoints("A")
ggsave(a,filename = "plots/dot_breakpoints_A.pdf", height = 6, width = 8)
g=make_dot_plot_for_breakpoints("G")
ggsave(g,filename = "plots/dot_breakpoints_G.pdf", height = 6, width = 8)
c=make_dot_plot_for_breakpoints("C")
ggsave(c,filename = "plots/dot_breakpoints_C.pdf", height = 6, width = 8)
t=make_dot_plot_for_breakpoints("T")
ggsave(t,filename = "plots/dot_breakpoints_T.pdf", height = 6, width = 8)
combining_dot_plot=ggarrange(a,g,c,t, ncol = 2, nrow = 2)
ggsave(combining_dot_plot,filename = "plots/combining_dot_breakpoints.pdf",height = 6, width = 8)
#example command line
#R --vanilla --slave --args sample_configure.csv 0 < /Users/chenhong1/DIRNA-StatsViz/DIRNA_vis_multiple_samples.R



