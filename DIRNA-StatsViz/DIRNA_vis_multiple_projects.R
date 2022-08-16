library(ggplot2)
library(ggpubr)
parameter=commandArgs(trailingOnly = TRUE)
projectconfigure=parameter[1]
projectconfigure=read.csv(projectconfigure,header=TRUE,sep= ",")
number_of_project=length(projectconfigure$Projectname)
make_dot_plot_for_breakpoints = function(Nucleotide)
{
  project1=read.table(paste(projectconfigure[1,"Projectname"],"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
  row.names(project1)=project1$Nucleotide
  project1_data=as.data.frame(t(project1[Nucleotide,2:5]))
  project1_data$Type=projectconfigure[1,1]
  project1_data$Position=row.names(project1_data)
  
  project2=read.table(paste(projectconfigure[2,"Projectname"],"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
  row.names(project2)=project2$Nucleotide
  project2_data=as.data.frame(t(project2[Nucleotide,2:5]))
  project2_data$Type=projectconfigure[2,"Projectname"]
  project2_data$Position=row.names(project2_data)
  master=rbind(project1_data,project2_data)
  for(projectname in projectconfigure[2:number_of_project,"Projectname"])
  {
    project=read.table(paste(projectname,"/breakpoints/breakpoints_observed_percentage.csv",sep=""),header=TRUE,sep=",")
    row.names(project)=project$Nucleotide
    project_data=as.data.frame(t(project[Nucleotide,2:5]))
    project_data$Type=projectname
    project_data$Position=row.names(project_data)
    master=rbind(master,project_data)
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

#R --vanilla --slave --args project_configure.csv < /Users/chenhong1/DIRNA-StatsViz/DIRNA_vis_multiple_projects.R