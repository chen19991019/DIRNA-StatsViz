The product DIRNA-StatsViz is designed accroding the composition of data and the type of analysis, like single samples analysis, multiple samples analysis and multiple projects analysis. The sacipts can be divided into four parts, such as the script for expected data of DI-RNAs breakpoints, single-sample scripts, multi-sample scripts and multi-project scripts.

# DIRNA_breakpoints_expected.py

It is the script for expected data of DI-RNAs breakpoints.

It is always called by single-sample scripts. It is responsible for assisting single-sample scripts in generating expected data like PB2_expected_ nucleotides.csv, breakpoints_expected_nucleotides.csv and breakpoints_expected_ percentage.csv of DI-RNAs for the sample. 

# DIRNA_stats_single_sample.py & DIRNA_vis_single_sample.R

They are single-sample scripts.and apllied to perform single-sample analysis.

## Input
You supposed to possess the analysed sample data folder and the .fasta files of eight segments of that sample.

## Usage

```
usage: DIRNA_stats_single_sample.py [-h] --sample SAMPLE 
                        [--noplots_for_singlesample]
                        [--project PROJECT]
                        [--cycle_index CYCLE_INDEX]

Analyse DI-RNAs in a single sample
optional arguments:
  -h, --help            show this help message and exit
  --sample SAMPLE       The name of the sample
  --noplots_for_singlesample
                        Use it only if you do not want to plot for the single
                        sample
  --project PROJECT     The name of the project
  --cycle_index CYCLE_INDEX
                        The number of cycles to produce the expected data for
                        breakpoints.Use it only if you want to customize the
                        number of cycle, the default number of cycle is equal
                        to the total count of the dirna of the sample
  ```
  After coming to the sample folder, the scripts are run by one command:
  Example command line using example data:
  ``./DIRNA-StatsViz/DIRNA_stats_single_sample.py --sample illumina``
  
 Scripts generate plots by default, but you can control the script not to draw plots and only obtain the desired data tables by adding option ```--noplots_for_singlesample``` . 
 the number of times for cycling in breakpoints scripts is by default equal to the total number of DI-RNAs counts of the sample, but the you are allowed to customise the number of cycles to n by ```--cycle_index``` . <n> is a number customised by the you.

## Output
Four data tables folders and a plot folder are produced to sample folder by default.

The data folders are:
* dirna_length
* dirna_length_percentage
* remaining_length
* breakpoints

The pdf of plots in plots folder are:
* combining_breakpoints.pdf
* combining_frequency_barplot.pdf
* combining_heatmap.pdf
* count_scale_plot.pdf
* frequency_of_nucleotides.pdf
* 8 heatmap_<segemnt>.pdf
* length_frequency_plot.pdf
* length_percentage_frequency_plot_.pdf
* remaining_length_frequency_plot_.pdf

## Warning
There will be an inevitable warning when drawing the WebLogo plot.For example:
  
Warning message:
`guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> = "none")` instead. 

# DIRNA_stats_multiple_samples.py & DIRNA_vis_multiple_samples.R
  
They are multi-sample scripts.
  
For analysing data of multiple samples, DIRNA_stats_multiple_samples.py and DIRNA_vis_multiple_samples.R can be used.

## Input
To analyse the DI-RNAs data of multiple samples, besides the data folder of each sample and .fasta FATSA files for eight segments refered by these samples, the user also needs to edit a sample_configure.csv(shown above) file for these samples and sort these files and folders to a new working folder finally.

## Usage

```
usage: DIRNA_stats_multiple_samples.py [-h] --dirna_count_threshold DIRNA_COUNT_THRESHOLD --venn_upset_threshold VENN_UPSET_THRESHOLD [--noplots_for_singlesample] [--noplots_for_multiplesamples]
                                       [--project PROJECT] [--cycle_index CYCLE_INDEX]

Analyse DI-RNAs between multiple samples of one project

optional arguments:
  -h, --help            show this help message and exit
  --dirna_count_threshold DIRNA_COUNT_THRESHOLD
                        The threshold used for filtering dirna by their counts
  --venn_upset_threshold VENN_UPSET_THRESHOLD
                        The threshold used to filtering by their counts when you plot venn and upset for multiple samples
  --noplots_for_singlesample
                        Use it only if you do not want to plot for the single sample
  --noplots_for_multiplesamples
                        Use it only if you do not want to plot for multiple samples
  --project PROJECT     The name of the project
  --cycle_index CYCLE_INDEX
                        The number of cycles to produce the expected data for breakpoints.Use it only if you want to customize the number of cycle, the default number of cycle is equal to the total
                        count of the dirna of the sample
  ```
  After coming to the working folder, the scripts are run by one command.
  Example command line using example data:
  ``python3 ./DIRNA-StatsViz/DIRNA_stats_multiple_samples.py --dirna_count_threshold 3 --venn_upset_threshold 3``
  
 Scripts both generate tables and plots for multiple samples and call single-sample scripts to analyse and plot for each single sample by default. You can add option ```--noplots_for_singlesample```and option ```--noplots_for_multiplesamples``` to prevent them generatring plots. 
 the number of times for cycling in breakpoints scripts is by default equal to the total number of DI-RNAs counts of the sample, but the you are allowed to customise the number of cycles to n by ```--cycle_index``` . <n> is a number customised by the you.

## Output
Four processing data tables folders, two useful tables and a plot folder are produced to working folder by default.

The data folders are
* combined_deletions
* filtered_deletions
* sorted_final_deletions
* breakpoints

Two useful tables are:
* Overview.csv
* Summary.csv

The pdf of plots in plots folder are:
* combining_line_for_DiRNA_type.pdf
* combining_line_for_DiRNA_count.pdf
* combining_venn.pdf
* 8 venn_<segemnt>.pdf
* 8 upset_<segemnt>.pdf
* dot_breakpoints_<nucleotide>.pdf
* combining_dot_breakpoints.pdf
* 8 line_count_<segemnt>.pdf
* summary_type.pdf
* combining_summary.pdf
* summary_count.pdf

## Warning
There will be an inevitable warning if sample sample has no qualifying DI-RNAs after further screening when drawing the Venn and UpSet plot.
  For example:
"!!!Warning:There is one sample where the count of all DiRNA of the NA segment is less than the threshold, so the upset plot does not apply to segment NA"
  
and 
  
geom_path: Each group consists of only one observation. Do you need to adjust
the group aesthetic? 

# DIRNA_stats_multiple_projects.py & DIRNA_vis_multiple_projects.R
  
They are multi-project scripts.
  
DIRNA_stats_multiple_projects.py and DIRNA_vis_multiple_projects.R were the scripts for anslysing DI-RNAs of multiple projects.

## Input
Before running the scripts, users are supposed to prepare a project_configure.csv (shown above) file for analysed projects firstly and put it, the .fasta files of the standard reference sequences and all project data folders together in a new working folder finally. Additionally, the sample_configure.csv file for multiple samples also needs to create in per project folder, as mentioned above.

## Usage

```
usage: DIRNA_stats_multiple_projects.py [-h] --dirna_count_threshold
                                        DIRNA_COUNT_THRESHOLD
                                        --venn_upset_threshold
                                        VENN_UPSET_THRESHOLD
                                        [--noplots_for_singlesample]
                                        [--noplots_for_multiplesamples]
                                        [--noplots_for_multipleprojects]
                                        [--cycle_index CYCLE_INDEX]

Analyse DI-RNAs between different projects

optional arguments:
  -h, --help            show this help message and exit
  --dirna_count_threshold DIRNA_COUNT_THRESHOLD
                        The threshold used for filtering dirna by their counts
  --venn_upset_threshold VENN_UPSET_THRESHOLD
                        The threshold used to filtering by their counts when
                        you plot venn and upset for multiple samples
  --noplots_for_singlesample
                        Use it only if you do not want to plot for the single
                        sample
  --noplots_for_multiplesamples
                        Use it only if you do not want to plot for multiple
                        samples
  --noplots_for_multipleprojects
                        Use it only if you do not want to plot for multiple
                        projects
  --cycle_index CYCLE_INDEX
                        The number of cycles to produce the expected data for
                        breakpoints.Use it only if you want to customize the
                        number of cycle, the default number of cycle is equal
                        to the total count of the dirna of the sample
  ```
  After coming to the working folder, the scripts are run by one command line.
  Example command line using example data:
  ``python3 ./DIRNA-StatsViz/DIRNA_stats_multiple_projects.py --dirna_count_threshold 3 --venn_upset_threshold 3``
  
 Scripts would call scripts for single sample and for multiple samples to be run, so they would generate plots for single sample, multiple samples and multiple projects by default.You can add option ```--noplots_for_singlesample```, option ```--noplots_for_multiplesamples``` and option ```--noplots_for_multipleprojects```to prevent them generatring plots. 
 the number of times for cycling in breakpoints scripts is by default equal to the total number of DI-RNAs counts of the sample, but the you are allowed to customise the number of cycles to n by ```--cycle_index``` . <n> is a number customised by the you.

## Output
Three processing data tables folders, a summary table and a plot folder are produced to working folder by default.

The data folders are:
* combined_deletions
* filtered_deletions
* sorted_final_deletions

A summary table is: Summary.csv

The pdf of plot in plots folder is: dot_breakpoints.pdf



## Warning
There will be an inevitable warning printed by single-sample scripts and multi-sample scripts.
