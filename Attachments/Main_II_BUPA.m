% 4th project of Teleriabilitazione course - part II.
% Zijlstra method - GaitRite initial contacts (ICs)/stride durations (SDs)
% comparison.
% Author: Burçak Carlo Pasqua.
% Contact: burcak.pasqua@gmail.com

clearvars
close all
clc 

% adding to the MATLAB search path the folder containing the handcrafted
% functions needed (in '02_functions' in the current project folder)
path_to_add = pwd;
path_to_add = strcat(path_to_add,filesep,'02_functions');
addpath(path_to_add);  % if old search path is wanted to be restored
                       % without the necessity to end the current Matlab
                       % session, just take the old Matlab session when
                       % updating it and restore it at the end with another
                       % call to 'addpath'

%% population choice and data loading

% pop choice
inpok = false;
while ~inpok
    pop = input("please choose the population which onto work:\n" + ...
        " - enter 'e' to work on healty elderly subjects\n" + ...
        " - enter 's' to work on stroke subjects\n\n choice: ","s");
    if pop == "e"  % healty elderly subjects
        pop = "Elderly";
        inpok = true;
    elseif pop == "s"  % stroke subjects 
        pop = "Stroke";
        inpok = true;
    else
        disp('please enter a valid choice')
    end
end

% data loading
fsep = filesep;
load ("."+fsep+"03_intermediate results & results"+fsep+"IC-"+pop+".mat")  % ICs instant estimates from MIMU data
load ("."+fsep+"01_data"+fsep+"GaitRite-"+pop+".mat")  % ICs instants from GaitRite gold standard data

%% ICs and SDs of GS matching estimates and descriptive statistics on their errors 

n_sub = length(GaitRite); 
descriptive_stat = struct;
for n = 1:n_sub
    
  % calculation of parameters of interest for current subject on matching 
  % detected gait events (GEs)
  [IC_MIMU,n_missed,n_extra,SD_MIMU,IC_GS,SD_GS] = IC_ON_GAITRITE( ...
      pop_parameters,GaitRite,n);
  
  % storage of current subject's parameters (descriptive statistics) in 
  % the global population structure
  descriptive_stat(n).SD_MIMU = SD_MIMU';
  descriptive_stat(n).SD_GS = SD_GS';
  IC_err = IC_MIMU-IC_GS; IC_err = IC_err(~isnan(IC_err)); % errors (in seconds) on ICs estimates
  descriptive_stat(n).IC_average_err = mean(IC_err);
  descriptive_stat(n).IC_std_err = std(IC_err);
  descriptive_stat(n).IC_MAE = mean(abs(IC_err));
  descriptive_stat(n).missed_percentage = n_missed/length(IC_GS)*100;
  descriptive_stat(n).extra_percentage = n_extra/length(IC_GS)*100;
  SD_err = SD_MIMU-SD_GS; SD_err = SD_err(~isnan(SD_err)); % errors (in seconds) on SDs estimates
  descriptive_stat(n).SD_average_err = mean(SD_err);
  descriptive_stat(n).SD_std_err = std(SD_err);
  descriptive_stat(n).SD_MAE = mean(abs(SD_err));
  
  % stride duration Shapiro-Wilk normality test & histogram
  H = SWTEST(descriptive_stat(n).SD_MIMU,0.05);  % normality test
  prc_5_95 = prctile(descriptive_stat(n).SD_MIMU,[5 95]);  % 5-th and 95-th percentiles
  if H == 1
    info = "distribution is normal according to Shapiro-Wilk " + ...
        "normality test";
  else
      info = "distribution is NOT normal according to Shapiro-Wilk " + ...
          "normality test";
      H = SWTEST(descriptive_stat(n).SD_MIMU(descriptive_stat(n).SD_MIMU ...
          >prc_5_95(1) & descriptive_stat(n).SD_MIMU<prc_5_95(2)),0.05);  % normality test after removal of data not included between 5-th and 95-th percentiles
      if H == 1
          info = info+", but it becomes normal after the removal of " + ...
              "data lower than 5-th percentile and greater than 95-th " + ...
              "percentile";
      else
          info = info+", and still remains the same after the " + ...
              "removal of data lower than 5-th percentile and " + ...
              "greater than 95-th percentile";
      end
  end
  figure,sgtitle("subject "+num2str(n)+" of "+pop+" population - " + ...
      "stride durations histogram")  % distribution representation
  ax = axes('OuterPosition',[.05 .05 .9 .75]);
  annotation("textbox",[.15 .85 .8 .05],'String',info,'EdgeColor','none')
  histogram(ax,descriptive_stat(n).SD_MIMU,round((max ...
      (descriptive_stat(n).SD_MIMU)-min(descriptive_stat(n).SD_MIMU))/ ...
      power(std(descriptive_stat(n).SD_MIMU,"omitnan"),2))), 
  xline(prc_5_95,'--',{'5-th percentile','95-th percentile'}, ...
      'LabelHorizontalAlignment','center')
  xlabel('stride duration (s)'),ylabel('occurences (#)')
 
end

% summary of descriptive statistics of errors on the population and, if
% evaluating errors on the stroke population, improvement evaluation with
% respect to results shown at slide 58 of project 4 provided outline
IC_AverageMAE = mean([descriptive_stat.IC_MAE]);  % average IC MAE on pop
SD_AveragePercentageMAE = mean([descriptive_stat.SD_MAE]./ ...  % average percentage SD MAE on pop
    arrayfun(@(x) mean(x.SD_GS),descriptive_stat))*100;  
fprintf(['\n\nSummary of descriptive statistics of errors for the ' ...
    'implemented method on %s population:\n\t- average MAE on ICs ' ...
    'detection: %f\n\t- average percentage MAE on SDs: %f%s\n'], ...
    pop,IC_AverageMAE,SD_AveragePercentageMAE,"%")
if pop == "Elderly"
    por = [0.010 3.9 1.1]; % provided outline results [MAE %missed %extra]
else
    por = [0.120 21 3.6]; % provided outline results [MAE %missed %extra]
end
if por(1)-IC_AverageMAE >= 0
    s1 = "-";
else
    s1 = "+";
end
missed_comp = por(2)-mean([descriptive_stat.missed_percentage]);  % comparison on average percentage of missed events
if  missed_comp >= 0
    s2 = "-";
else
    s2 = "+";
end
extra_comp = por(3)-mean([descriptive_stat.extra_percentage]);  % comparison onaverage percentage of extra events
if extra_comp >= 0
    s3 = "-";
else
    s3 = "+";
end
fprintf(['\nWith respect to results presented in the given outline,' ...
        'the implemented method shows:\n\t- IC MAE: %s%s\n\t- missed ' ...
        'ICs: %s%s\n\t- extra ICs: %s%s\n'],s1,string(abs((por(1)- ...
        IC_AverageMAE)/por(1))*100)+"%",s2,string(abs(missed_comp))+ ...
        "%",s3,string(abs(extra_comp))+"%")

%% MIMU estimates vs GaitRite goldstandard Bland-Altman plot 

% whole population SDs
pop_SD_MIMU = [descriptive_stat.SD_MIMU];
pop_SD_GS = [descriptive_stat.SD_GS];

% keeping only SD of correctly found ICs
pop_SD_GS = pop_SD_GS(~isnan(pop_SD_MIMU));
pop_SD_MIMU = pop_SD_MIMU(~isnan(pop_SD_MIMU));

% Shapiro-Wilk normality test on all the differences of SDs values of the 
% population, this is important to know if we can assume mean ± 1.96*SD 
% represents the 95% confidence interval or not
pop_SD_dbm = pop_SD_MIMU-pop_SD_GS;  % stride duration paired differences 
                                     % between MIMU based and GS methods
H = SWTEST(pop_SD_dbm,0.05);  % normality test
if H == 1
    division = false;  % flag for division of data (see plot below)
    info = "distribution is normal according to Shapiro-Wilk test, " + ...
        "therefore mean ± 1.96*st.dev. respresents the 95% " + ...
        "confidence interval";
    m = mean(pop_SD_dbm);
    stdev = std(pop_SD_dbm);
    RI = [m-1.96*stdev m m+1.96*stdev]; % reference interval
else
    info = "distribution is NOT normal according to Shapiro-Wilk " + ...
        "normality test";
    prc_5_95 = prctile(pop_SD_dbm,[5 95]);
    H = SWTEST(pop_SD_dbm(pop_SD_dbm>prc_5_95(1)&pop_SD_dbm<prc_5_95(2)),0.05);  % normality test after removal of data not included between 5-th and 95-th percentiles
      if H == 1
          division = true;  % flag for division of data (see plot below)
          info = info+", but it becomes normal considering only " + ...
              "data included between the 5-th and the 95-th " + ...
              "percentiles, therefore even if the data are all " + ...
              "plotted the mean ± 1.96*st.dev. are referred only to " + ...
              "these last so that it represents their 95% " + ...
              "confidence interval";
          pop_SD_dbm_in595 = pop_SD_dbm(pop_SD_dbm>prc_5_95(1) & pop_SD_dbm<prc_5_95(2));
          pop_SD_dbm_out595 = pop_SD_dbm(pop_SD_dbm<prc_5_95(1) | pop_SD_dbm>prc_5_95(2));
          m = mean(pop_SD_dbm_in595);
          stdev = std(pop_SD_dbm_in595);
          RI = [m-1.96*stdev m m+1.96*stdev]; % reference interval
      else
          % (for stroke population) in case the distribution of differences
          % between the MIMU estimates and GS still is not normal, this is 
          % probably due to the highly impaired subjects 1 and 7 of the
          % stroke population on which the method works with poorer 
          % performances than on the other stroke population subjects, 
          % anyway removal of these patients from further considerations is
          % not done as the size of the sample is too small to consider
          % them as outliers (while notice that removing data outer from 
          % 5-th and 95-th percentiles is not a same entity issue as these 
          % could be only some of the data of subjects for which the method 
          % puntually did not performed well during the trial, eventually 
          % due to many reasons)
          division = false;  % flag for division of data (see plot below)
          info = info+", and still remains the same even " + ...
              "considering only data included between the 5-th and " + ...
              "the 95-th percentiles, therefore 2.5-th, 50-th and " + ...
              "97.5-th percentiles are represented rather than mean " + ...
              "± 1.96*st.dev.";
          RI = prctile(pop_SD_dbm,[5 50 95]);  % reference interval
      end
end
RI_lab = [string(RI(1)), string(RI(2)) string(RI(3))];  % labels

% Bland_Altman plot
figure,sgtitle('Bland-Altman plot - MIMU estimates SD vs GS SD ')
ax = axes('OuterPosition',[.05 .05 .9 .75]);
annotation("textbox",[.15 .85 .8 .05],'String',info,'EdgeColor','none')
if division
    plot(ax,pop_SD_GS(pop_SD_dbm>prc_5_95(1) & pop_SD_dbm<prc_5_95(2)), ...
        pop_SD_dbm_in595,'k.',pop_SD_GS(pop_SD_dbm<prc_5_95(1) | ...
        pop_SD_dbm>prc_5_95(2)),pop_SD_dbm_out595,'c.')
    legend('included in 5th-95th percentiles range',['excluded from ' ...
        '5th-95th percentiles range'])
else
    plot(ax,pop_SD_GS,pop_SD_dbm,'k.')
end
yline([RI(1) RI(3)],'k--',{RI_lab(1) RI_lab(3)}, ...
    'LabelVerticalAlignment','middle'),
yline(RI(2),'k',RI_lab(2),'LabelVerticalAlignment','middle'),
xlabel('stride duration GS (s)'),
ylabel('stride duration MIMU - stride duration GS'),
xlim([.5 5]),ylim([-1 1])

%% error descriptive statiscs saving
save(string(pwd)+fsep+"03_intermediate results & results"+fsep+ ...
    "ErrorsDescripriteStat-"+pop+".mat","descriptive_stat", ...
    "IC_AverageMAE","SD_AveragePercentageMAE")
