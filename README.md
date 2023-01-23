# Gait-temporal-parameters-estimation-from-accelerometry-data
The purpose was to implement a method proposed in literature to detect initial contacts of the walking motor task on forward direction
accelerometry signals acquired from a single accelerometer at pelvis level, with particular interest in thorny signals acquired from
hemiplegic subjects. Furthermore an innovative algorithm remarkably outperforming the aforementioned in terms of accuracy has been
proposed. Read REPORT.pdf for more detailed technical description of the project.

This project is licensed under the terms of the MIT license.

## Attached files
In ./Attachments/ subfolder (please read the report -REPORT.pdf- before for a better understanding of attachments' description):
* Main_I_BUPA.m: main implementing the initial contact instants and stride duration estimation;
* Main_II_BUPA.m: main implementing the performance analysis of the found estimates with respect to used gold standard (GaitRite electronic walkway); 
* 01_data/: subfolder containing the .mat format files of the time series signals under analysis;
* 02_functions/ZILSTRA_METHOD_BUPA.m: function implementing the method proposed in literature;
* 02_functions/icburc.m: function implementing the novel method;
* 02_functions/IC_ON_GAITRITE.p: function implementing calculation of parameters of interest (initial contacts and stride durations) on matching detected gait events;
* 02_functions/SWTEST.p: function implementing Shapiro-Wilk normality test;
* 03_intermediate results & results/: subfolder containing the .mat format files with initial contact estimates and related errors descrition.
