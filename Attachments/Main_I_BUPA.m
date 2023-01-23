% 4th project of Teleriabilitazione course - part I.
% Temporal gait parameters estimation from signals from a single MIMU on 
% pelvis.
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

%% population choice

fsep = filesep;  % file separator on the currently running on platform
dfrp = "." + fsep + "01_data" + fsep;  % relative path to folder with data
inpok = false;
while ~inpok
    pop = input("please choose the population which onto work:\n" + ...
        " - enter 'e' to work on healty elderly subjects\n" + ...
        " - enter 's' to work on stroke subjects\n\n choice: ","s");
    if pop == "e"
        ids = 1:10;  % healty elderly subjects IDs 1 -> 10 
        inpok = true;
        pop = "Elderly";  % population
    elseif pop == "s"
        ids = 21:29;  % stroke subjects IDs 21 -> 29
        inpok = true;
        pop = "Stroke";  % population
    else
        disp('please enter a valid choice')
    end
end

%% ICs and strides duration estimates for the population of interest

pop_parameters = struct;  % initialization of global population structure
cs = 0;  % counter on current subject of the population
for sub = ids  % loop on each subject of the population of interest
    
    % load of current subject's data
    if sub < 10
        fn = "S00";
    else
        fn = "S0";
    end
    fn = fn + num2str(sub) + "-dataMIMU.mat";
    load(dfrp+fn)
    cs = cs+1;
    
    % data unpacking (signals are already referred to standard anatomical 
    % RS )
    a_AP = MIMU.LowerBack.acc(:,3);  % anterior-posterior acceleration (with 
                                     % reference to anatomic orientation) 
    fs = MIMU.Fs;  % sampling frequency

    % ICs estimates (s) with peak detection method prposed by Zilstra et 
    % al. in 2003, notice that acceleration will not be de-normalized by g 
    % as the method performance are the same up to a scale factor
    a_AP = a_AP-mean(a_AP);  % mean value removal (otherwise proper 
                             % zero-crossing detection may fail)... for
                             % real this is no longer a problem in the new
                             % implemented ethod in icburc, BUT shall be
                             % done as otherwise if offset is no null a
                             % high power contribute can be carried by the
                             % constant component frequency performing an
                             % errate LP envelope extraction filtering 
                             % (even if anyway offset removal in embedded 
                             % in the function...)
    % IC_t = ZIJLSTRA_METHOD_BUPA(a_AP,[],fs,[0.1],[2],[20],0,"time");
    if pop == "Elderly"
        IC_t = icburc(a_AP,fs,0.1,17,TimeConversion="on");
    else
        IC_t = icburc(a_AP,fs,0.2,19,Type="asymm",TimeConversion="on");
    end
    
    % stride duration calculation assuming an IC for each gait cycle has 
    % been found (otherwise alternantion bewtween right and left foot 
    % strikes -which is necessarily present by definition of walking motor 
    % task- cannot be exploited to easily find stride durations starting 
    % from ICs detected)                                   
    if mod(length(IC_t),2)==1
        IC_t = [IC_t;NaN];  %#ok<AGROW> (message supression as it is not the case)
    end
    str_d = reshape(IC_t',2,[])'; str_d = diff(str_d); str_d = str_d(~isnan(str_d));
  
    % storage of current subject's parameters in the global population structure
    pop_parameters(cs).IC = IC_t;
    pop_parameters(cs).Stride_Duration = str_d;

end 

%% estimates data saving

save(string(pwd)+fsep+"03_intermediate results & results"+fsep+"IC-"+...
    pop+".mat","pop_parameters")

%% population check plot 

figure, hold on 
for sub = ids
    SD = pop_parameters(sub-ids(1)+1).Stride_Duration;
    plot(sub,SD,'.b'), 
end 
title('Stride Duration distribution across ' + pop + ' subjects'),
xlabel('# subject'), ylabel('Time (s)'),xlim([ids(1)-1 ids(end)+1])