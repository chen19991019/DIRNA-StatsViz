import argparse
import csv
import os
import subprocess
import random
import pandas as pd

# create the argument parser
parser = argparse.ArgumentParser(description='Analyse DI-RNAs in a single sample')
parser.add_argument('--sample',required=True, help='The name of the sample')
parser.add_argument('--noplots_for_singlesample',required=False, action='store_true',help='Use it only if you do not want to plot for the single sample')
parser.add_argument('--project',help='The name of the project')
parser.add_argument('--cycle_index',help='The number of cycles to produce the expected data for breakpoints.Use it only if you want to customize the number of cycle, '
                                         'the default number of cycle is equal to the total count of the dirna of the sample')
args = parser.parse_args()

if args.sample:
    print("Start data analysis for a single sample!")
    if args.project:
        #When performing analysis between projects, the path to each sample folder will change, so the file path will need to be changed.
        master_path =  args.project + "/" + args.sample + "/"
    else:
        master_path = args.sample + "/"
    newpath = master_path + "dirna_length"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    newpath = master_path + "dirna_length_percentage"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    newpath = master_path + "remaining_length"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    newpath = master_path + "breakpoints"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    #generate dirna_length.csv
    FILE = {"{}_PB2_dels_map.txt".format(args.sample): "PB2_dirna_length.csv",
            "{}_PB1_dels_map.txt".format(args.sample): "PB1_dirna_length.csv",
            "{}_PA_dels_map.txt".format(args.sample): "PA_dirna_length.csv",
            "{}_HA_dels_map.txt".format(args.sample): "HA_dirna_length.csv",
            "{}_NP_dels_map.txt".format(args.sample): "NP_dirna_length.csv",
            "{}_NA_dels_map.txt".format(args.sample): "NA_dirna_length.csv",
            "{}_M_dels_map.txt".format(args.sample): "M_dirna_length.csv",
            "{}_NS_dels_map.txt".format(args.sample): "NS_dirna_length.csv"}
    for key in FILE:
        try:
            file_path = master_path + key
            with open(file_path, 'r') as textFile:
                outfilename="{}dirna_length/{}".format(master_path,FILE[key])
                with open(outfilename, 'w', newline="") as out:
                    writer = csv.writer(out)
                    next(textFile)
                    n = 0
                    n1 = 0
                    n2 = 0
                    n3 = 0
                    n4 = 0
                    n5 = 0
                    n6 = 0
                    n7 = 0
                    n8 = 0
                    for line in textFile:
                        l = line.split('\t')
                        n += int(l[1])#get the total number of
                        if int(l[4]) >= 1 and int(l[4]) < 50:  # DI-RNA length
                            n1 += int(l[1])
                        elif int(l[4]) >= 50 and int(l[4]) < 100:
                            n2 += int(l[1])
                        elif int(l[4]) >= 100 and int(l[4]) < 250:
                            n3 += int(l[1])
                        elif int(l[4]) >= 250 and int(l[4]) < 500:
                            n4 += int(l[1])
                        elif int(l[4]) >= 500 and int(l[4]) < 1000:
                            n5 += int(l[1])
                        elif int(l[4]) >= 1000 and int(l[4]) < 1500:
                            n6 += int(l[1])
                        elif int(l[4]) >= 1500 and int(l[4]) < 2000:
                            n7 += int(l[1])
                        elif int(l[4]) >= 2000:
                            n8 += int(l[1])
                    try:
                        LINE = [["Length", "Count", "Percentage"],
                                ["1<=length<50", n1, round(n1 / n, 4)],
                                ["50<=length<100", n2, round(n2 / n, 4)],
                                ["100<=length<250", n3, round(n3 / n, 4)],
                                ["250<=length<500", n4, round(n4 / n, 4)],
                                ["500<=length<1000", n5, round(n5 / n, 4)],
                                ["1000<=length<1500", n6, round(n6 / n, 4)],
                                ["1500<=length<2000", n7, round(n7 / n, 4)],
                                ["2000<=length", n8, round(n8 / n, 4)]]
                    except ZeroDivisionError:
                        continue
                    for line in LINE:
                        writer.writerow(line)
        except FileNotFoundError:
            print("{} not found while generate percentage file. Please check and try again!".format(file_path))
            raise SystemExit(1)
    #generate dirna_length_percentage.csv
    FILE_LENGTH = {"{}_HA_dels_map.txt".format(args.sample): "HA_dirna_length_percentage.csv",
                   "{}_PA_dels_map.txt".format(args.sample): "PA_dirna_length_percentage.csv",
                   "{}_PB1_dels_map.txt".format(args.sample): "PB1_dirna_length_percentage.csv",
                   "{}_PB2_dels_map.txt".format(args.sample): "PB2_dirna_length_percentage.csv",
                   "{}_M_dels_map.txt".format(args.sample): "M_dirna_length_percentage.csv",
                   "{}_NA_dels_map.txt".format(args.sample): "NA_dirna_length_percentage.csv",
                   "{}_NP_dels_map.txt".format(args.sample): "NP_dirna_length_percentage.csv",
                   "{}_NS_dels_map.txt".format(args.sample): "NS_dirna_length_percentage.csv"}
    for key in FILE_LENGTH:
        try:
            file_path = master_path + key
            with open(file_path, 'r') as textFile:
                outfilename = "{}dirna_length_percentage/{}".format(master_path,FILE_LENGTH[key])
                with open(outfilename, 'w', newline="") as out:
                    writer = csv.writer(out)
                    next(textFile)
                    p = 0
                    p1 = 0
                    p2 = 0
                    p3 = 0
                    p4 = 0
                    p5 = 0
                    p6 = 0
                    p7 = 0
                    p8 = 0
                    p9 = 0
                    p10 = 0
                    for line in textFile:
                        l = line.split('\t')
                        p += int(l[1]) #get the total number of
                        if float(l[5]) >= 0 and float(l[5]) < 0.2:
                            p1 += int(l[1])
                        elif float(l[5]) >= 0.2 and float(l[5]) < 0.4:  # length
                            p2 += int(l[1])
                        elif float(l[5]) >= 0.4 and float(l[5]) < 0.5:
                            p3 += int(l[1])
                        elif float(l[5]) >= 0.5 and float(l[5]) < 0.6:
                            p4 += int(l[1])
                        elif float(l[5]) >= 0.6 and float(l[5]) < 0.7:
                            p5 += int(l[1])
                        elif float(l[5]) >= 0.7 and float(l[5]) < 0.8:
                            p6 += int(l[1])
                        elif float(l[5]) >= 0.8 and float(l[5]) < 0.9:
                            p7 += int(l[1])
                        elif float(l[5]) >= 0.9 and float(l[5]) <= 1:
                            p8 += int(l[1])
                    try:
                        LINE_LENGTH = [["Length_Percentage", "Count", "Percentage"],
                                       ["0<=length%<20", p1, round(p1 / p, 4)],
                                       ["20<=length%<40", p2, round(p2 / p, 4)],
                                       ["40<=length%<50", p3, round(p3 / p, 4)],
                                       ["50<=length%<60", p4, round(p4 / p, 4)],
                                       ["60<=length%<70", p5, round(p5 / p, 4)],
                                       ["70<=length%<80", p6, round(p6 / p, 4)],
                                       ["80<=length%<90", p7, round(p7 / p, 4)],
                                       ["90<=length%<=100", p8, round(p8 / p, 4)]]
                    except ZeroDivisionError:
                        continue
                    for line_length in LINE_LENGTH:
                        writer.writerow(line_length)
        except FileNotFoundError:
            print("{} not found while generate length_percentage file. Please check and try again!".format(key))
            raise SystemExit(1)
    #generate remaining_length.csv
    try:
        file_path = master_path + "{}_coverage.txt".format(args.sample)
        with open(file_path, 'r') as coverageFile:
            next(coverageFile)
            coverage = {}
            for line in coverageFile:
                line = line.replace('\t', ',')
                l = line.split(',')
                coverage[l[0]] = l[1]
        FILE_REMAIN_LENGTH = {"{}_HA_dels_map.txt".format(args.sample): {coverage["HA"]: "HA_remaining_length.csv"},
                            "{}_PA_dels_map.txt".format(args.sample): {coverage["PA"]: "PA_remaining_length.csv"},
                            "{}_PB1_dels_map.txt".format(args.sample): {coverage["PB1"]: "PB1_remaining_length.csv"},
                            "{}_PB2_dels_map.txt".format(args.sample): {coverage["PB2"]: "PB2_remaining_length.csv"},
                            "{}_M_dels_map.txt".format(args.sample): {coverage["M"]: "M_remaining_length.csv"},
                            "{}_NA_dels_map.txt".format(args.sample): {coverage["NA"]: "NA_remaining_length.csv"},
                            "{}_NP_dels_map.txt".format(args.sample): {coverage["NP"]: "NP_remaining_length.csv"},
                            "{}_NS_dels_map.txt".format(args.sample): {coverage["NS"]: "NS_remaining_length.csv"}}
        for key in FILE_REMAIN_LENGTH:
            try:
                file_path = master_path + key
                with open(file_path, 'r') as textFile:
                    for key2 in FILE_REMAIN_LENGTH[key]:
                        outfilename = "{}remaining_length/{}".format(master_path,FILE_REMAIN_LENGTH[key][key2])
                        with open(outfilename, 'w', newline="") as out:
                            writer = csv.writer(out)
                            next(textFile)
                            L = 0
                            l1 = 0
                            l2 = 0
                            l3 = 0
                            l4 = 0
                            l5 = 0
                            l6 = 0
                            l7 = 0
                            l8 = 0
                            l9 = 0
                            l10 = 0
                            for line1 in textFile:
                                line1 = line1.replace('\t', ',')
                                l = line1.split(',')
                                remainingcoordinate = l[0].split('-')
                                L += int(l[1])
                                #remaining length= coverage - DelEnd + DelStar
                                remaining = int(key2) - int(remainingcoordinate[1]) + int(remainingcoordinate[0])
                                if remaining >= 1 and remaining < 50:  # filter whether there is a dirna
                                    l1 += int(l[1])
                                elif remaining >= 50 and remaining < 100:  # length
                                    l2 += int(l[1])
                                elif remaining >= 100 and remaining < 250:
                                    l3 += int(l[1])
                                elif remaining >= 250 and remaining < 500:
                                    l4 += int(l[1])
                                elif remaining >= 500 and remaining < 1000:
                                    l5 += int(l[1])
                                elif remaining >= 1000 and remaining < 1500:
                                    l6 += int(l[1])
                                elif remaining >= 1500 and remaining < 2000:
                                    l7 += int(l[1])
                                elif remaining >= 2000:
                                    l8 += int(l[1])
                            try:
                                LINE = [["Remaining_Length", "Count", "Percentage"],
                                        ["1<=length<50", l1, round(l1 / L, 4)],
                                        ["50<=length<100", l2, round(l2 / L, 4)],
                                        ["100<=length<250", l3, round(l3 / L, 4)],
                                        ["250<=length<500", l4, round(l4 / L, 4)],
                                        ["500<=length<1000", l5, round(l5 / L, 4)],
                                        ["1000<=length<1500", l6, round(l6 / L, 4)],
                                        ["1500<=length<2000", l7, round(l7 / L, 4)],
                                        ["2000<=length", l8, round(l8 / L, 4)]]
                            except ZeroDivisionError:
                                continue
                            for line in LINE:
                                writer.writerow(line)
            except FileNotFoundError:
                print("{} cannot be open while generate remaining_length file!".format(key))
                raise SystemExit(1)
    except FileNotFoundError:
        print("{}_coverage.txt cannot be open while generate remaining_length file!".format(args.sample))
        raise SystemExit(1)

#generate breakpoints_observed data for each sample
if args.project:
    files = os.listdir(args.project+'/'+args.sample)
else:
    files = os.listdir(args.sample)
datas_up = []
datas_down = []
sumcount=0
for file in files:
    if '_dels_map.txt' in file:
        filepath=master_path+file
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
                datas_up.append(data_up)#put all up_dirna sequence in to datas_up
                datas_down.append(data_down)#put all down_dirna sequence in to datas_down
        sumcount += count
    else:
        pass

A_up = [0] * 2
G_up = [0] * 2
C_up = [0] * 2
T_up = [0] * 2
A_down = [0] * 2
G_down = [0] * 2
C_down = [0] * 2
T_down = [0] * 2
for data in datas_up:
    for i in range(0, 2):
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
A=A_up + A_down
A=[round(a/sumcount,4) for a in A] # keep 4 decimal places
G=G_up + G_down
G=[round(g/sumcount,4) for g in G]
C=C_up + C_down
C=[round(c/sumcount,4) for c in C]
T=T_up + T_down
T=[round(t/sumcount,4) for t in T]
head = ["Nucleotide", "up_2", "up_1", "down_1", "down_2"]
with open("{}breakpoints/breakpoints_observed_percentage.csv".format(master_path), 'w', newline="") as breakpoint:
    writer = csv.writer(breakpoint)
    writer.writerow(head)
    writer.writerow(["A"] + A)
    writer.writerow(["G"] + G)
    writer.writerow(["C"] + C)
    writer.writerow(["T"] + T)

#build pipeline to generate breakpoints_expected data
script_path = os.path.dirname(__file__)
if args.project:
    if args.cycle_index:
        commandline='python3 {}/DIRNA_breakpoints_expected.py --master_path {} --cycle_index {} --project {}'.format(
            script_path,master_path,args.cycle_index,args.project)
    else:#The default is to generate the expected data based on the loop with the total number of dirna counts of all segments of the sample.
        commandline = 'python3 {}/DIRNA_breakpoints_expected.py --master_path {} --cycle_index {} --project {}'.format(
            script_path, master_path, sumcount, args.project)
else:
    if args.cycle_index:
        commandline = 'python3 {}/DIRNA_breakpoints_expected.py --master_path {} --cycle_index {}'.format(
            script_path, master_path, args.cycle_index)
    else:
        commandline = 'python3 {}/DIRNA_breakpoints_expected.py --master_path {} --cycle_index {}'.format(
            script_path,master_path,sumcount)
run_expected = subprocess.run(commandline,shell=True)


# build the pipeline to run R
if args.noplots_for_singlesample==False:
    print("Start visualization for a single sample using R!")
    newpath = master_path + "plots"
    if not os.path.exists(newpath):
        os.makedirs(newpath)
    script_path = os.path.dirname(__file__)
    if args.project:
        workdirectory =  args.project
    else:
        workdirectory = ''
    Rcomandline = "R --vanilla --slave --args {} {} < {}/DIRNA_vis_single_sample.R".format(args.sample, workdirectory,script_path)
    runR=subprocess.run(Rcomandline, shell=True)
    print("All analyses for a single sample finished!")
else:
    print("Data analysis for a single sample finished!")

#some example command lines
#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_single_sample.py  --sample minion --cycle_index 1000 --noplots_for_singlesample
#python3 /Users/chenhong1/DIRNA-StatsViz/DIRNA_stats_single_sample.py  --sample illumina --project Hutchinson
