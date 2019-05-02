% Filename: EMG_processing_tutorial.m
% Author:   Samuel Acuña
% Date:     11 May 2018
% Description:
% Guides you through a typical EMG processing routine.
%
% EMG is more about timing of activation, than the amplitude.
% 
% Make sure you have the accompanying data file:
% EMG_processing_tutorial_DATA.mat  or
% EMG_processing_tutorial_DATA2.mat
%
% Directions: run each section at a time to see steps (usually alt+enter, or command+enter)

%% STEP 1: LOAD EMG DATA
clear; close all; clc;
%load('EMG_processing_tutorial_DATA.mat');
load('EMG_processing_tutorial_DATA2.mat');
% DATA: this is data from the medial gastrocnemius when walking on a treadmill
% DATA2: tibialis anterior, same trial

% Loads matlab structure: trial_data
%   trial_data.emg         : raw EMG data
%   trial_data.label       : muscle where EMG data was collected
%   trial_data.freq        : sampling frequency
%   trial_data.time        : time of each sample of data
%   trial_data.heelStrikes : instances of heel strikes between steps

return
%% STEP 1: EXAMINE RAW EMG DATA
figure(1);
h1 = plot(trial_data.time, trial_data.emg,'color',[  0    0.4470    0.7410]);
title(['EMG: ' trial_data.label]);
xlabel('time (sec)'); xlim([0 60]);
ylabel('EMG (volts)');
legend([h1],'raw EMG');

%% STEP 2: VISUALIZE HEEL STRIKES
figure(1); hold on;
scale = max(trial_data.emg);
hs_x = [trial_data.heelStrikes, trial_data.heelStrikes]';
hs_y = [zeros(length(trial_data.heelStrikes),1)-(scale/2), zeros(length(trial_data.heelStrikes),1)+(scale/2)]';
h2 = plot(hs_x,hs_y,'r');
hold off;
legend([h1 h2(1)],{'raw EMG', 'heel strikes'});

%% STEP 3: BANDPASS FILTER EMG
% the bandpass gets rid of any noise that isnt part of EMG data. Most EMG
% power is between 5-500 Hz
% low cutoff: want to be high enough to get rid of motion artifacts and drift. No higher than 10 hz, per ISEK guidelines
% high cutoff: no lower than 350 Hz, per ISEK guidelines
BP = [10 500]; % bandpass filter parameters in Hz [low cutoff, high cutoff]
[b_BP,a_BP]=butter(4,BP/(trial_data.freq/2)); % bandpass filter coefficients. (4th order butterworth)
emg_BP = filtfilt(b_BP,a_BP,trial_data.emg); % zero phase lag filter

figure(1); hold on;
h3 = plot(trial_data.time,emg_BP,'color',[0.6875    0.7656    0.8672]);
hold off;
legend([h1 h2(1) h3],{'raw EMG', 'heel strikes', 'bandpass EMG'});
%% STEP 4: FULL-WAVE RECTIFY EMG
% just the absolute value
emg_ABS = abs(emg_BP);

figure(1); hold on;
h4 = plot(trial_data.time,emg_ABS,'color',[0.6758    0.8438    0.8984]);
hold off;
legend([h1 h2(1) h3 h4],{'raw EMG', 'heel strikes', 'bandpass EMG', 'rectify EMG'});

%% STEP 5: LINEAR ENVELOPE
% low pass filter
LP = 10; % low pass filter for linear envelope, in Hz
[b_LP,a_LP]=butter(4,LP/(trial_data.freq/2),'low'); % linear envelope filter
emg_ENV = filtfilt(b_LP,a_LP,emg_ABS);

figure(1); hold on;
h5 = plot(trial_data.time,emg_ENV,'color',[0  0 0],'LineWidth',2);
hold off;
legend([h1 h2(1) h3 h4 h5],{'raw EMG', 'heel strikes', 'bandpass EMG', 'rectify EMG', 'envelope EMG'});


%% STEP 6: NORMALIZE AMPLITUDE
% scaled to the max value (could also do RMS, or max voluntary contraction)
emg_NORM = emg_ENV/max(emg_ENV);

figure(2);
h6 = plot(trial_data.time,emg_NORM,'color',[ 0    0.4470    0.7410],'LineWidth',2);
title(['EMG: ' trial_data.label]);
xlabel('time (sec)'); xlim([0 60]);
ylabel('normalized EMG'); ylim([0 1.5]);
legend([h6],'normalized linear envelope of EMG');

%% STEP 7: TIME NORMALIZE BY A STRIDE

% divides the EMG signals insto strides
npts = 101; % points per gait cycle
nStrides = length(trial_data.heelStrikes)-1; % number of complete gait cycles
emg_Strides = zeros(101,nStrides); %preallocate
for j = 1:nStrides
    j1 = find(trial_data.time>trial_data.heelStrikes(j),1); % get index of time at first heel strike
    j2 = find(trial_data.time>trial_data.heelStrikes(j+1),1); % get index of time at second heel strike
    emg_Strides(:,j) = normcycle(emg_NORM(j1:j2),npts); % time normalize
end

% FIND AVERAGE EMG over a stride
emg_AVG = mean(emg_Strides,2); %average EMG of each stride
emg_STD = std(emg_Strides')'; % standard deviation


% PLOT
figure(3); subplot(2,1,1); 
shadedErrorBar([0:100]',emg_AVG,emg_STD);
title(['EMG: ' trial_data.label]);
xlabel('Gait Cycle (0-100%)'); xlim([0 100]);
ylabel('normalized EMG'); ylim([0 1.5]);
legend('EMG (AVG ± STD)')

figure(3); subplot(2,1,2); hold on;
for i = 1:nStrides
        plot([0:100]',emg_Strides(:,i));
end
hold off;
title(['EMG for every stride']);
xlabel('Gait Cycle (0-100%)'); xlim([0 100]);
ylabel('normalized EMG'); ylim([0 1.5]);



function yf = normcycle(y,n,x)
% yf = normcycle(y,n,x)
% Convert a signal y to n even-spaced data points over a cycle
% Often used for presentation of gait data, default for n is 101 points
% can specify an indpendent variable x (optional)
if ~exist('n','var')
    n=101;
end
[nr,nc]=size(y);
if nc==1 && nr>1
    ny=1;
    nx=nr;
elseif nr==1 && nc>1
    y=y';
    ny=1;
    nx=nc;
elseif nr>1 && nc>1
    ny=nc;
    nx=nr;
else
    disp('normcycle does not work on a scalar value');
    yf=[];
    return
end
if ~exist('x','var')
    x=[0:(nx-1)]/(nx-1);
else
    nx=length(x);
    x=(x-x(1))/(x(end)-x(1));
end
kk=[0:(n-1)]/(n-1);
yf=interp1(x,y,kk,'*pchip');

end
