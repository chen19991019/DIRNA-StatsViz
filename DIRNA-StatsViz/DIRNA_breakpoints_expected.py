#%%
import argparse
import random
import csv
import pandas as pd

# create the argument parser
parser = argparse.ArgumentParser(description='Generate the expected data for breakpoints')
parser.add_argument('--master_path', help='The master path of the working folder')
parser.add_argument('--project', help='The name of the project')
parser.add_argument('--cycle_index',required=True,help='The number of cycles to produce the expected data for breakpoints')
args = parser.parse_args()

datas=[]
data_up_e=[]
data_down_e=[]
files=[]
SEGMENT = ["PB2", "PB1", "PA", "HA", "NP", "NA", "M", "NS"]
for se in SEGMENT:
    for t in range(0, int(args.cycle_index)):
        if args.project:
            fastafile='{}/{}.fasta'.format(args.project,se)
        else:
            fastafile = '{}.fasta'.format(se)
        with open(fastafile, 'r') as textFile:
            next(textFile)
            n = 0
            for line in textFile:
                seqlist = list(line)
                RefLength = len(line)
        RandomNumber1 = random.randint(20, RefLength - 20 - 25)
        RandomNumber2 = random.randint(RandomNumber1, RefLength - 20)
        up = []
        down = []
        for i in range(2, 0, -1):
            up.append(seqlist[RandomNumber1 - i])
        for i in range(1, 3, 1):
            down.append(seqlist[RandomNumber2 + i])
        data_up_e.append(up)
        data_down_e.append(down)
    for i in range(0, len(data_up_e)):
        data = data_up_e[i] + data_down_e[i]
        datas.append(data)
    head = ["up_2", "up_1", "down_1", "down_2"]
    if args.master_path:
        expected_data_filename="{}breakpoints/{}_expected_nucleotides.csv".format(args.master_path,se)
    else:
        expected_data_filename = "breakpoints/{}_expected_nucleotides.csv".format(se)
    with open(expected_data_filename, 'w', newline="") as breakpoint:
        writer = csv.writer(breakpoint)
        writer.writerow(head)
        writer.writerows(datas)
    f = pd.read_csv(expected_data_filename)
    files.append(f)
File = pd.concat(files)
if args.master_path:
    expectedfile="{}breakpoints/breakpoints_expected_nucleotides.csv".format(args.master_path)
else:
    expectedfile="breakpoints/breakpoints_expected_nucleotides.csv"
File.to_csv(expectedfile, index=0)
up_1=[]
up_2=[]
down_1=[]
down_2=[]

with open(expectedfile, 'r') as dataFile:
    next(dataFile)
    for line in dataFile:
        line=line.strip("\n")
        line=line.split(',')
        up_2.append(line[0])
        up_1.append(line[1])
        down_1.append(line[2])
        down_2.append(line[3])
if args.master_path:
    expected_percentage_file='{}breakpoints/breakpoints_expected_percentage.csv'.format(args.master_path)
else:
    expected_percentage_file='breakpoints/breakpoints_expected_percentage.csv'
with open(expected_percentage_file, 'w', newline="") as finalout:
    writer = csv.writer(finalout)
    head = ["Nucleotide","up_2", "up_1", "down_1", "down_2"]
    writer.writerow(head)
    for Nu in ('A', 'G', 'C', 'T'):
        row=[Nu,round(up_2.count(Nu)/len(up_2),4),round(up_1.count(Nu)/len(up_1),4),
             round(down_1.count(Nu)/len(down_1),4),round(down_2.count(Nu) / len(down_2),4)]
        writer.writerow(row)


