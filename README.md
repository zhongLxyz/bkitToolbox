The bkitToolbox is a data analysis toolbox written in MATLAB designed for the Biomotion Lab to analyze dataset output from PsychBench and BMLkit. It is responsible for most of the analysis need in the study of the perceptual bias in the perception of biological motion. It can visualize and model the data, generate new data to simulate the characteristics of a desired participant among other things.

Here is the documentation for the functions included in the toolbox.

- bkitRead(name)

this function reads in data files from all participants (sessions), the output is in
 the forms of a cell array, with each cell containing
 a table of  data for one participant. input the csv file name or the folder name to process the data in that folder.

- bkitFigure
This function plots the relevant graphs for bmlkit data.

arguments: 
how to use it: the first argument has to be the dataset read from bkitRead. The other arguments are optional,
 including 

 (1) the name of the plot you want. 'FTVelevation', 'VFAelevation', 'VFAazimuth' for the 3 kinds of plots.

   FTVelevation: for each session, plot the percentage of FTV responses at each elevation angles.

   VFAelevation: for each session, plot the percentage of VFA responses at each elevation angles.

   VFAazimuth: for each session, plot the percentage of VFA responses at each azimuth angles.

 (2) specification of what kind of responses is coded in the response column of the dataset.
'ftv', 'cw', 'vfa' - it specifies what the kind of responses is coded in the response column.

- plot3DAvg
This script does the following things:

1. creating 3D bar graph of bias size for each participant in each condition (3*3 bar graph). including an average plot for all participants.

2. printing a table containing the average bias size across all participants for each conditions, along with their SD

3. print a table with t-test results comparing whether there is a significant difference between
%    the size of (the same) biases in different conditions.

- bkitGenerate(parameter, design, allBalancedVar, ntrial)

This function generates simulation of a experiment session using parameter values
 fitted to the GLM from bkitModel, along with the designs of that experiment from bkitSummary. It will calculate the simulated response of a 'participant' from 
the given values.

- bkitModel(allData, type, varargin)

This function finds the model parameters from glm for the specified datasets, it basically finds what variable setting in a walker has an effect on the perception of the participants, and how big the effect is.
 for each individual dataset, it outputs a set of model parameters.

- bkitSummary(allData)

This function reads in the output of bkitRead (the output from all participants) and outputs a summary
 for the dataset. The summary includes the design of the experiment. The design is what variables are in the dataset and what levels they have. A dataset passed to it
 can have multiple designs in it. It will also find what variables are balanced together
in the dataset

- bkitSequence
This script finds whether there is any serial dependency effect in the input dataset.
It will conduct the test in all 3 variables (cw/ccw, ftv/fa, vfa/vfb)
It conducts a t-test for the data of each session, to determine if
repeating answers happen more often than chance would predict.

- NbackCorrelationCW
This function reads in the dataset from a CW/CCW condition and finds the
correlation between each trial and its previous 1~N trials (see below for
more information)
This is one measure of the serial dependency effect.

- NbackCorrelationFTV
This function reads in the dataset from a FTV/FA condition and finds the
correlation between each trial and its previous 1~N trials (see below for
more information)
This is one measure of the serial dependency effect.

- NbackCorrelationVFA
This function reads in the dataset from a viewing perspective condition and finds the
correlation between each trial and its previous 1~N trials (see below for
more information)
This is one measure of the serial dependency effect. 

- NbackFormatted
This script is created to aid in creating the table for a summary. 
It generates the formatted result of 1- to N-back serial dependency analysis of
the dataset, using the format of the summary table, the result T can be copied into the
summary directly.

- findCCW
This function finds the trials in which the responses are CCW.input: angularVelocity, azimuth, elevation and response columns in the dataset. The
response column can be coded in 0 and 1, or -1 and 1, Angular Velocity
can be coded in 45 and -45, or 1 and -1.
output: a column vector in which FTV=1 and FA=0 for corresponding trials.