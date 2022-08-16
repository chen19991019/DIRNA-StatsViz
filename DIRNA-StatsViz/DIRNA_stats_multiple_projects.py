import argparse
import os
import subprocess
import pandas as pd
import csv
import shutil
#%%
# create the argument parser
parser = argparse.ArgumentParser(description='Analyse DI-RNAs between different projects')
#parser.add_argument('--projects', required=True,nargs='+',help='The name of different projects')
parser.add_argument('--dirna_count_threshold', required=True, type = int, help='The threshold used for filtering dirna by their counts')
parser.add_argument('--venn_upset_threshold', type = int,required=True, help='The threshold used to filtering by their counts when you plot venn and upset for multiple samples')
parser.add_argument('--noplots_for_singlesample',required=False, action='store_true',help='Use it only if you do not want to plot for the single sample')
parser.add_argument('--noplots_for_multiplesamples',required=False, action='store_true',help='Use it only if you do not want to plot for multiple samples')
parser.add_argument('--noplots_for_multipleprojects',required=False, action='store_true',help='Use it only if you do not want to plot for multiple projects')
parser.add_argument('--cycle_index',help='The number of cycles to produce the expected data for breakpoints.Use it only if you want to customize the number of cycle, '
                                         'the default number of cycle is equal to the total count of the dirna of the sample')
args = parser.parse_args()

projects=[]
with open("project_configure.csv", "r") as project_configureFile:
        next(project_configureFile)
        for line in project_configureFile:
            line = line.strip('\n')
            line = line.split(',')
            projects.append(line[0])

print("Start data analysis for multiple samples, a single sample and finally multiple projects!")
for pro in projects:
    multiplescript_path = os.path.dirname(__file__)
    if args.noplots_for_multiplesamples == False and args.noplots_for_singlesample == False:
        if args.cycle_index:
            commandline ='python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --cycle_index {}'.format(
            multiplescript_path,pro,args.dirna_count_threshold,args.venn_upset_threshold,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {}'.format(
                multiplescript_path,pro,args.dirna_count_threshold,args.venn_upset_threshold)
    elif args.noplots_for_singlesample == False and args.noplots_for_multiplesamples == True:
        if args.cycle_index:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --cycle_index {} --noplots_for_multiplesamples'.format(
            multiplescript_path,pro,args.dirna_count_threshold,args.venn_upset_threshold,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --noplots_for_multiplesamples'.format(
                multiplescript_path,pro,args.dirna_count_threshold,args.venn_upset_threshold)
    elif args.noplots_for_multiplesamples == False and args.noplots_for_singlesample == True:
        if args.cycle_index:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --cycle_index {} --noplots_for_singlesample'.format(
            multiplescript_path, pro, args.dirna_count_threshold, args.venn_upset_threshold,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --noplots_for_singlesample'.format(
                multiplescript_path, pro, args.dirna_count_threshold, args.venn_upset_threshold)
    elif args.noplots_for_multiplesamples == True and args.noplots_for_singlesample == True:
        if args.cycle_index:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --cycle_index {} --noplots_for_singlesample --noplots_for_multiplesamples'.format(
            multiplescript_path, pro, args.dirna_count_threshold, args.venn_upset_threshold,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_multiple_samples.py --project {} --dirna_count_threshold {} --venn_upset_threshold {} --noplots_for_singlesample --noplots_for_multiplesamples'.format(
                multiplescript_path, pro, args.dirna_count_threshold, args.venn_upset_threshold)
    out = subprocess.run(commandline, shell=True)
    
print("Start data analysis for multiple projects!")
SEGMENT = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
newpath = "prepare"
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = "combined_deletions"
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = "filtered_deletions"
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = "final_deletions"
if not os.path.exists(newpath):
    os.makedirs(newpath)
newpath = "sorted_final_deletions"

if not os.path.exists(newpath):
    os.makedirs(newpath)
for seg in SEGMENT:
    FILEDATA = []
    files=[]
    for pro in projects:
        columnname = []
        SAMPLELIST = []
        filename = "{}/sorted_final_deletions/{}_sorted_final_deletions.csv".format(pro,seg)
        file = pd.read_csv(filename)
        for name in file:
            columnname.append(name)
        LINE = []
        n = 0
        with open(filename, "r") as starfile:
            next(starfile)
            for line in starfile:
                line = line.strip('\n')
                line = line.split(',')
                samplelists = '{} : '.format(pro)
                for a in range(len(line) - 1, 5, -1):
                    if int(line[a]) != 0:
                        samplelist = '{}={}|'.format(columnname[a], line[a]) #only after I add - behind them ,san I seperate the list不然，那些样品和数量都挤在一起
                        samplelists=samplelists+samplelist
                samplelists=samplelists.strip("|")
                line = [line[0], line[1], line[2], line[3]] + [samplelists]
                LINE.append(line)
            csvname = "prepare/{}_{}.csv".format(pro,seg)
            csv_head = ["Segment", "DelStart-DelEnd", "Count", "Number_of_Sample", "List_of_Samples"]
            with open(csvname, 'w', newline="") as combinedcsv:
                writer = csv.writer(combinedcsv)
                writer.writerow(csv_head)
                writer.writerows(LINE)
        FILE = "prepare/{}_{}.csv".format(pro, seg)
        f = pd.read_csv(FILE)
        files.append(f)
    File = pd.concat(files)
    combinedname="combined_deletions/{}_combined_deletions.csv".format(seg)
    File.to_csv(combinedname, index=0)
    filter = pd.read_csv(combinedname)
    filter.drop_duplicates(subset=['DelStart-DelEnd'], inplace=True)
    filtercsvname = "filtered_deletions/{}_filtered_deletions.csv".format(seg)
    filter.to_csv(filtercsvname, index=0)

for seg1 in SEGMENT:
    DATA=[]
    filtercsvname = "filtered_deletions/{}_filtered_deletions.csv".format(seg1)
    with open(filtercsvname, "r") as filteredFile:
        L = []
        next(filteredFile)
        for LINE in filteredFile:
            count = 0
            number_of_sampples = 0
            list_of_samples = str()
            LINE = LINE.split(',')
            for pro in projects:
                FILE = "prepare/{}_{}.csv".format(pro, seg1)
                with open(FILE, "r") as textFile:
                    lines = []
                    next(textFile)
                    for line in textFile:
                        line = line.strip('\n')
                        line = line.split(',')
                        if LINE[1] == line[1]:
                            count += int(line[2])
                            number_of_sampples += int(line[3])
                            list_of_samples=list_of_samples+str(line[4])+' ; '
            list_of_samples=list_of_samples.strip(" ; ")
            data = [LINE[0], LINE[1], number_of_sampples, count]+[list_of_samples]
            DATA.append(data)

    csvname = "final_deletions/{}_final_deletions.csv".format(seg1)
    csv_head = ["Segment", "DelStart-DelEnd", "Number_of_Sample","Count", "List_of_Samples"]
    with open(csvname, 'w', newline="") as combinedcsv:
        writer = csv.writer(combinedcsv)
        writer.writerow(csv_head)
        writer.writerows(DATA)
    sortfile = pd.read_csv(csvname)
    sortfile = sortfile.sort_values(["Number_of_Sample", "Count"],
                                    ascending=[False, False])  # filter by multiple columns
    sortfile.to_csv("sorted_final_deletions/{}_sorted_final_deletions.csv".format(seg1), index=0)

#Summary table
files=[]
for seg in SEGMENT:
    sorted_final_deletions_name = "sorted_final_deletions/{}_sorted_final_deletions.csv".format(seg)
    f=pd.read_csv(sorted_final_deletions_name)
    files.append(f)
File = pd.concat(files)
File = File.sort_values(["Number_of_Sample", "Count"],ascending=[False, False])
File.to_csv("Summary.csv", index=0)


#breakpoints_observed for each project
#Combine all the dirna in a project into one whole.
for pro in projects:
    samples = []
    datas_up = []
    datas_down = []
    sumcount = 0
    master_path = pro + "/"
    sample_configureFile_path = master_path + "sample_configure.csv"
    with open(sample_configureFile_path, "r") as sample_configureFile:
        next(sample_configureFile)
        for line in sample_configureFile:
            l = line.split(',')
            samples.append(l[0])
    for sample in samples:
        files = os.listdir(master_path+sample)
        for file in files:
            if '_dels_map.txt' in file:
                filepath = master_path + sample + '/' + file
                count = 0
                with open(filepath, "r") as textFile:
                    next(textFile)
                    for line in textFile:
                        line = line.strip("\n")
                        l = line.split('\t')
                        count += int(l[1])
                        data_up = list(l[2])
                        data_up.append(l[1])
                        data_down = list(l[3])
                        data_down.append(l[1])
                        datas_up.append(data_up)
                        datas_down.append(data_down)
                sumcount += count
            else:
                pass
    A_up = [0] * 2  # 5个的时候，换成5
    G_up = [0] * 2
    C_up = [0] * 2
    T_up = [0] * 2
    A_down = [0] * 2
    G_down = [0] * 2
    C_down = [0] * 2
    T_down = [0] * 2
    for data in datas_up:
        for i in range(0, 2):  # 5个的时候换成5
            if data[i] == 'A':
                A_up[i] += int(data[2])
            elif data[i] == 'G':
                G_up[i] += int(data[2])
            elif data[i] == 'C':
                C_up[i] += int(data[2])
            elif data[i] == 'T':
                T_up[i] += int(data[2])
    for data in datas_down:
        for i in range(0, 2):
            if data[i] == 'A':
                A_down[i] += int(data[2])
            elif data[i] == 'G':
                G_down[i] += int(data[2])
            elif data[i] == 'C':
                C_down[i] += int(data[2])
            elif data[i] == 'T':
                T_down[i] += int(data[2])
    try:
        A = A_up + A_down
        A = [round(a / sumcount, 4) for a in A]
        G = G_up + G_down
        G = [round(g / sumcount, 4) for g in G]
        C = C_up + C_down
        C = [round(c / sumcount, 4) for c in C]
        T = T_up + T_down
        T = [round(t / sumcount, 4) for t in T]
    except ZeroDivisionError:
        continue

    head = ["Nucleotide", "up_2", "up_1", "down_1", "down_2"]
    breakpointsname = "{}breakpoints/breakpoints_observed_percentage.csv".format(master_path)
    with open(breakpointsname, 'w', newline="") as breakpoint:
        writer = csv.writer(breakpoint)
        writer.writerow(head)
        writer.writerow(["A"] + A)
        writer.writerow(["G"] + G)
        writer.writerow(["C"] + C)
        writer.writerow(["T"] + T)
print("Data analysis for multiple projects finished!")

#Whether you want to plot for the multiple projects
if args.noplots_for_multipleprojects==False:
    print("Start visualization for multiple projects using R!")
    newpath = "plots"
    if not os.path.exists(newpath):
       os.makedirs(newpath)
    script_path = os.path.dirname(__file__)
    Rcomandline = "R --vanilla --slave --args project_configure.csv < {}/DIRNA_vis_multiple_projects.R".format(script_path)
    out = subprocess.run(Rcomandline,shell=True)
    print("All analyses for multiple projects finished!")
else:
    print("All analyses for multiple projects finished!")
#delete useless folder
shutil.rmtree("final_deletions")
shutil.rmtree("prepare")


#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_multiple_projects.py --projects 1012_2-3-6 1012_2-3 1013_2-3 --dirna_count_threshold 0 --venn_upset_threshold 0 --cycle_index 10000 --noplots_for_singlesample  --noplots_for_multiplesamples --noplots_for_multipleprojects
