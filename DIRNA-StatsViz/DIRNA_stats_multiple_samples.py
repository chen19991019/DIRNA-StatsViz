import argparse
import os
import subprocess
import pandas as pd
import csv
import shutil#use to delete useless folder
# create the argument parser
parser = argparse.ArgumentParser(description='Analyse DI-RNAs between multiple samples of one project')
#parser.add_argument('--sample_configureFile', required=True, help='The configure file containing sample names and some specific information')
parser.add_argument('--dirna_count_threshold', required=True, type = int, help='The threshold used for filtering dirna by their counts')
parser.add_argument('--venn_upset_threshold', type = int,required=True, help='The threshold used to filtering by their counts when you plot venn and upset for multiple samples')
parser.add_argument('--noplots_for_singlesample',required=False, action='store_true',help='Use it only if you do not want to plot for the single sample')
parser.add_argument('--noplots_for_multiplesamples',required=False, action='store_true',help='Use it only if you do not want to plot for multiple samples')
parser.add_argument('--project', help='The name of the project')
parser.add_argument('--cycle_index',help='The number of cycles to produce the expected data for breakpoints.Use it only if you want to customize the number of cycle, '
                                         'the default number of cycle is equal to the total count of the dirna of the sample')
args = parser.parse_args()

samples=[]

print("Start data analysis for multiple samples!")
if args.project:
    # When performing analysis between projects, the path to each sample folder will change, so the file path will need to be changed.
    master_path =  args.project + "/"
    newpath1 = master_path + "combined_deletions"
    newpath2 = master_path + "filtered_deletions"
    newpath3 = master_path + "final_deletions"
    newpath4 = master_path + "sorted_final_deletions"
    newpath5 = master_path + "breakpoints"
    sample_configureFile_path = master_path + "sample_configure.csv"
else:
    newpath1 = "combined_deletions"
    newpath2 = "filtered_deletions"
    newpath3 = "final_deletions"
    newpath4 = "sorted_final_deletions"
    newpath5 = "breakpoints"
    sample_configureFile_path = "sample_configure.csv"
if not os.path.exists(newpath1):
    os.makedirs(newpath1)
if not os.path.exists(newpath2):
    os.makedirs(newpath2)
if not os.path.exists(newpath3):
    os.makedirs(newpath3)
if not os.path.exists(newpath4):
    os.makedirs(newpath4)
if not os.path.exists(newpath5):
    os.makedirs(newpath5)
with open(sample_configureFile_path, "r") as sample_configureFile:
    next(sample_configureFile)
    for line in sample_configureFile:
        line=line.strip('\n')
        l = line.split(',')
        samples.append(l[0])
SEGMENT = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
#Merge the same fragment of .dels_map.txt file and remove the duplicated dirna coordinates, only leaving unique coordinates
for seg in SEGMENT:
    FILEDATA = []
    for a in range(0, len(samples), 1):
        if args.project:
            FILE = ["{}{}/{}_{}_dels_map.txt".format(master_path,samples[a], samples[a], seg)]
            FILEDATA.append(FILE)
        else:
            FILE = ["{}/{}_{}_dels_map.txt".format(samples[a], samples[a], seg)]
            FILEDATA.append(FILE)
    datas = []
    for file in FILEDATA:
        file_path = file[0]
        with open(file_path, "r") as textFile:
            next(textFile)
            for line in textFile:
                line = line.replace('\t', ',')
                l = line.split(',')
                datas.append(l)
    if args.project:
        csvname = "{}combined_deletions/{}_combined_deletions.csv".format(master_path,seg)
    else:
        csvname = "combined_deletions/{}_combined_deletions.csv".format(seg)
    csv_head = ["DelStart-DelEnd", "Count", "5-2bp", "3-2bp", "Length", "Length%"]
    with open(csvname, 'w', newline="") as combinedcsv:
        writer = csv.writer(combinedcsv)
        writer.writerow(csv_head)
        writer.writerows(datas)
    #remove the duplicated dirna coordinates, only leaving unique coordinates
    filter = pd.read_csv(csvname)
    filter.drop_duplicates(subset=['DelStart-DelEnd'], inplace=True)
    if args.project:
        filtercsvname = "{}filtered_deletions/{}_filtered_deletions.csv".format(master_path,seg)
    else:
        filtercsvname = "filtered_deletions/{}_filtered_deletions.csv".format(seg)
    filter.to_csv(filtercsvname, index=0, encoding='gbk')
#Using the filtered unique coordinates as reference, loop through each file to find the coordinate that is same as them
#and get the related information of that coordinate
for seg1 in SEGMENT:
    if args.project:
        filtercsvname = "{}filtered_deletions/{}_filtered_deletions.csv".format(master_path,seg1)
    else:
        filtercsvname = "filtered_deletions/{}_filtered_deletions.csv".format(seg1)
    L = []
    with open(filtercsvname, "r") as filteredFile:
        next(filteredFile)
        for LINE in filteredFile:
            Samplecount = [0] * len(samples)  # prepare to record number of deletion possessed by each sample
            n = 0
            count = 0
            filter_count = 0
            LINE = LINE.split(',')
            for a1 in range(0, len(samples), 1):
                if args.project:
                    FILE = "{}{}/{}_{}_dels_map.txt".format(master_path, samples[a1], samples[a1], seg1)
                else:
                    FILE = "{}/{}_{}_dels_map.txt".format(samples[a1], samples[a1], seg1)
                with open(FILE, "r") as textFile:
                    next(textFile)
                    for line in textFile:
                        line = line.replace('\t', ',')
                        line = line.split(',')
                        if LINE[0] == line[0]:
                            n += 1
                            count += int(line[1])
                            Samplecount[a1] = int(line[1])
                            segment = seg1
                            if segment == "NA":  # !!!NA will show as blank in the sorted and overview table's column, so change to Na
                                segment = "Na"
                            if int(line[1]) > args.dirna_count_threshold:  #the count of the specific diRNA is greater than dirna_count_threshold;
                                # dirna_count_threshold could be 3 or 10 or any user defined value.
                                filter_count += 1
            average = round(count / n, 3)
            L1 = [segment, LINE[0], count, n, filter_count, average] + Samplecount
            L.append(L1)
    head1 = ["Segment", "DelStart-DelEnd", "Count", "Number_of_Sample",
                "Number_of_Sample(Count>{})".format(args.dirna_count_threshold),
                "Average_Count"]
    for a2 in range(0, len(samples), 1):
        h2 = "{}".format(samples[a2])
        head1.append(h2)
    if args.project:
        final_deletions_name = "{}final_deletions/{}_final_deletions.csv".format(master_path, seg1)
    else:
        final_deletions_name = "final_deletions/{}_final_deletions.csv".format(seg1)
    with open(final_deletions_name, 'w', newline="") as finalout1:
        writer = csv.writer(finalout1)
        writer.writerow(head1)
        writer.writerows(L)
    sortfile = pd.read_csv(final_deletions_name)
    sortfile = sortfile.sort_values(["Number_of_Sample", "Average_Count"],
                                    ascending=[False, False])  # filter by multiple columns
    if args.project:
        sorted_final_deletions_name = "{}sorted_final_deletions/{}_sorted_final_deletions.csv".format(master_path, seg1)
    else:
        sorted_final_deletions_name = "sorted_final_deletions/{}_sorted_final_deletions.csv".format(seg1)
    sortfile.to_csv(sorted_final_deletions_name, index=0)

#%%summary and overview table
files=[]
L = []
ssum_conut = [0] * len(samples)
sshared_conut = [0] * len(samples)
ssum_dirna = [0] * len(samples)
for seg in SEGMENT:
    if args.project:
        sorted_final_deletions_name = "{}sorted_final_deletions/{}_sorted_final_deletions.csv".format(master_path, seg)
    else:
        sorted_final_deletions_name = "sorted_final_deletions/{}_sorted_final_deletions.csv".format(seg)
    f=pd.read_csv(sorted_final_deletions_name)
    files.append(f)
    overalcount=0
    overal_dirna = 0
    share_count=0
    share_dirna = 0
    notshare_count=0
    notshare_dirna = 0
    share_count_filter = 0
    share_dirna_filter = 0
    with open(sorted_final_deletions_name, "r") as textFile:
        next(textFile)
        for line in textFile:
            l = line.split(',')
            overalcount += int(l[2])
            overal_dirna += 1
            if int(l[3]) == len(samples):
                share_count += int(l[2])
                share_dirna += 1
            else:
                notshare_count += int(l[2])
                notshare_dirna += 1
            if int(l[4]) == len(samples) and int(l[3]) == len(samples):  # those whose count >=dirna_count_threshold
                share_count_filter += int(l[2])
                share_dirna_filter += 1
            se = seg
            if se == "NA":
                se = "Na"
    for a1 in range(0, len(samples), 1):
        with open(sorted_final_deletions_name, "r") as textFile:
            next(textFile)
            samplesum_conut = 0
            samplesum_dirna = 0
            sampleshared_conut = 0
            scindex = int(a1 + 6)
            for line in textFile:
                l = line.split(',')
                samplesum_conut += int(l[scindex])#sample's deletion
                if int(l[scindex]) > 0:
                    samplesum_dirna += 1
                if int(l[3]) == len(samples):
                    sampleshared_conut += int(l[scindex])
        ssum_conut[a1] = samplesum_conut
        sshared_conut[a1] = sampleshared_conut
        ssum_dirna[a1] = samplesum_dirna
    L2 = [se, overalcount, overal_dirna, share_count, share_dirna, share_count_filter, share_dirna_filter, notshare_count, notshare_dirna]+ssum_conut+sshared_conut+ssum_dirna
    L.append(L2)

#Summary table
head3 = ["Segment","Sum_Counts", "Sum_Types", "Shared_Counts", "Shared_Types", "Shared_Counts(Count>{})".format(args.dirna_count_threshold), "Shared_Types(Count>{})".format(args.dirna_count_threshold), "Not_Shared_Counts", "Not_Shared_Types"]
for a2 in range(0, len(samples), 1):
    h2 = "{}(Sum_Counts)".format(samples[a2])
    head3.append(h2)
for a3 in range(0, len(samples), 1):
    h3 = "{}(Shared_Counts)".format(samples[a3])
    head3.append(h3)
for a4 in range(0, len(samples), 1):
    h4 = "{}(Sum_Types)".format(samples[a4])
    head3.append(h4)
if args.project:
    summary_name = "{}Summary.csv".format(master_path)
else:
    summary_name = "Summary.csv"
with open(summary_name, 'w', newline="") as finalout1:
    writer = csv.writer(finalout1)
    writer.writerow(head3)
    writer.writerows(L)

#overview table-merge the sorted_final_deletions file of eight segments
File = pd.concat(files)
if args.project:
    overview_name = "{}Overview.csv".format(master_path)
else:
    overview_name = "Overview.csv"
File.to_csv(overview_name, index=0)
samples_for_one=[]
with open(sample_configureFile_path, "r") as sample_configureFile:
    next(sample_configureFile)
    for line in sample_configureFile:
        l = line.split(',')
        samples_for_one.append(l[0])




# build the pipeline to analyse each sample
for s in samples_for_one:
    # Data analysis for a single sample completed and start visualization using R
    onescript_path = os.path.dirname(__file__)
    if args.project:
        if args.cycle_index:
            commandline = 'python3 {}/DIRNA_stats_single_sample.py --sample {} --project {} --cycle_index {}'.format(
                onescript_path,s,args.project,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_single_sample.py --sample {} --project {}'.format(
                onescript_path, s, args.project)
    else:
        if args.cycle_index:
            commandline = 'python3 {}/DIRNA_stats_single_sample.py --sample {} --cycle_index {}'.format(
                onescript_path,s,args.cycle_index)
        else:
            commandline = 'python3 {}/DIRNA_stats_single_sample.py --sample {}'.format(
                onescript_path, s)
    # Whether you want to plot for the multiple samples
    if args.noplots_for_singlesample == False:
        out = subprocess.run(commandline, shell=True)
    else:
        out = subprocess.run(commandline + ' --noplots_for_singlesample', shell=True)


#breakpoints_common_expected
#The same reference sequence is used for multiple samples,
# so,they share a expected data.100000 cycles are used to produce the expected data when performing multiple-sample analysis.
script_path = os.path.dirname(__file__)
if args.project:
    commandline = 'python3 {}/DIRNA_breakpoints_expected.py --master_path {} --cycle_index {} --project {}'.format(
        script_path, master_path, 100000, args.project)
else:
    commandline = 'python3 {}/DIRNA_breakpoints_expected.py --cycle_index {}'.format(
        script_path, 100000)
run_expected = subprocess.run(commandline,shell=True)

print("Data analysis for multiple samples finished!")


#Whether you want to plot for the multiple samples
if args.noplots_for_multiplesamples==False:
    print("Start visualization for multiple samples using R!")
    if args.project:
        newpath5 = args.project + "/" + "plots"
    else:
        newpath5 = "plots"
    if not os.path.exists(newpath5):
       os.makedirs(newpath5)
    script_path = os.path.dirname(__file__)
    if args.project:
        workdirectory = args.project
    else:
        workdirectory = ''
    Rcomandline = "R --vanilla --slave --args sample_configure.csv {} {} < ".format(args.venn_upset_threshold,workdirectory)
    out = subprocess.run(Rcomandline + script_path + '/DIRNA_vis_multiple_samples.R', shell=True)
    print("All analyses for multiple samples finished!")
else:
    print("All analyses for multiple samples finished!")

#delete useless folder
if args.project:
    shutil.rmtree(args.project+"/final_deletions")
else:
    shutil.rmtree("final_deletions")

#some example command lines
#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_multiple_samples.py --dirna_count_threshold 3 --venn_upset_threshold 0 --cycle_index 100 --project Project2
#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_multiple_samples.py --dirna_count_threshold 3 --venn_upset_threshold 3 --noplots_for_multiplesamples --noplots_for_singlesample
#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_multiple_samples.py --dirna_count_threshold 3 --venn_upset_threshold 3 --project Project2
