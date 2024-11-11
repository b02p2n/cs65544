java c
Digital Biosignal Processing
MATLAB Laboratory 2
The objective of this exercise is to familiarise yourselves with the concept of discrete time convolution and discrete time Fourier transform. by generating an EMG signal from experimental motor unit action potentials (MUAPs) and motor neuron discharge timings. EMG signals can be modelled as a convolutive mixture of a series of delta functions representing the discharge timings of motor neurons in the spinal cord. The impulse responses of this convolutive mixture are the action potentials of the muscle units (Figure below; see slides 16-21 of Lecture 2).

For this laboratory, you are provided with the discharge timings and MUAPs recorded from the biceps brachii muscle of a healthy individual during a contraction at constant torque. The MUAPs are stored in the file “MUAPs.mat”, the discharge timings (as series of discrete time delta functions) in the file “NeuralDrive.mat”, and the recorded torque in the file “Torque.mat”. The sampling frequency of the recordings is 2048 Hz.
Study the short Matlab script. provided below. Use the Matlab script. to analyse the discharge timings and MUAPs in the frequency domain (Fourier transform). Use the script. also for reconstructing the EMG signal and analysing its Fourier transform.
The impulse response of a moving average filter is provided in the script, together with the plot of the frequency response of the filter. Analyse the frequency response by varying the length of the filter. Finally, complete the script. with the computation of the EMG envelope by filtering the rectified EMG (check the effect of the filter length on the envelope estimate).
In your report, please provide the following:- Use the Fourier transform. of the discharge timings to estimate the average discharge rate of the first motor neuron and comment on the way you obtained the estimate. [30%]- A plot of torque and EMG envelopes for filter lengths of 1000, 5000, and 10000. Comment on the effect of filter duration in relation to the frequency response of the filter and to the estimated envelope. [70%]
PLEASE NOTICE: The report is limited to one A4 page, including all graphs and comments.

% Clear working space
clear all
close all
clc
% Load required signals
load('MUAPs.mat'); % Single motor unit action potentials (experimental)
load('NeuralDrive.mat'); % Discharge times of motor neurons (experimental)
load('Torque.mat'); % Experimental Torque
fsamp = 2048; % Sampling frequency of the recordings
%% PART 1: Reconstructing the EMG signal through convolution of discharge times by MUAPs
n_MUAPs = size(MUAPs,1); % Number of MUAPs
dur_MUAPs = size(MUAPs,2); % Duration of MUAPs
dur_MUAPseq = size(Real_firing(1,:),2); % Duration of the signal
time_ax=0:1/fsamp:(dur_MUAPseq-1)/fs% Time axis for the signal
% Plot MUAP trains
figure(1),
for jj = 1:n_MUAPs
conv_train = conv(Real_firing(jj,:),MUAPs(jj,:)); % Convolution (see slide 14-15 of Lecture 2)
MUAP_Train(jj,:) = co代 写Digital Biosignal Processing MATLAB Laboratory 2Matlab
代做程序编程语言nv_train(floor(dur_MUAPs/2)+1:end-floor(dur_MUAPs/2)); % Cut transitory portion
hold on, plot(time_ax,MUAP_Train(jj,:)/(max(MUAPs(:)) - min(MUAPs(:))) + (n_MUAPs - jj + 1)*1.25,'k');
end
title('Sequence of MUAPs for each motor unit')
xlabel('Time (s)')
ylabel ('MUAP trains')
MUAP_sel = 1; % Select one MUAP (1 to 15) for Fourier analysis
% Plot Fourier Transform. of the discharge timings
figure(2)
f_transf_Firings = fft(Real_firing(MUAP_sel,:));
freq_ax = [-pi+pi/dur_MUAPseq:2*pi/dur_MUAPseq:pi-pi/dur_MUAPseq];
plot(freq_ax,fftshift(abs(f_transf_Firings)));
xlabel('Discrete Angular Frequency')
title('Discrete time Fourier Transform. of Motor Neuron Discharge Sequence')
ylabel('Magnitude of Fourier Transform. (Arbitrary Units)')
% Plot Fourier Transform. of the MUAP
figure(3)
f_transf_MUAP = fft(MUAPs(MUAP_sel,:));
freq_ax = [-pi+pi/dur_MUAPs:2*pi/dur_MUAPs:pi-pi/dur_MUAPs];
plot(freq_ax,fftshift(abs(f_transf_MUAP)));
xlabel('Discrete Angular Frequency')
title('Discrete time Fourier Transform. of Motor Unit Action Potential')
ylabel('Magnitude of Fourier Transform. (Arbitrary Units)')
% Plot Fourier Transform. of the MUAP train
figure(4)
f_transf_MUAP_Train = fft(MUAP_Train(MUAP_sel,:));
freq_ax = [-pi+pi/dur_MUAPseq:2*pi/dur_MUAPseq:pi-pi/dur_MUAPseq];
plot(freq_ax,fftshift(abs(f_transf_MUAP_Train)));
xlabel('Discrete Angular Frequency')
title('Discrete time Fourier Transform. of Motor Unit Action Potential Train')
ylabel('Magnitude of Fourier Transform. (Arbitrary Units)')
% Obtaining the EMG signal by summing all MUAP trains
recoEMG = sum(MUAP_Train,1);
figure(5), plot(1/fsamp:1/fsamp:length(recoEMG)/fsamp,recoEMG);
title('Reconstructed EMG signal'), xlabel('Time (s)'); ylabel('EMG (Arbitrary Units)');
% Plot Fourier Transform. of the EMG signal
figure(6)
f_transf_EMG = fft(recoEMG(MUAP_sel,:));
freq_ax = [-pi+pi/dur_MUAPseq:2*pi/dur_MUAPseq:pi-pi/dur_MUAPseq];
plot(freq_ax,fftshift(abs(f_transf_EMG)));
xlabel('Discrete Angular Frequency')
title('Discrete time Fourier Transform. of the EMG Signal')
ylabel('Magnitude of Fourier Transform. (Arbitrary Units)')
%% PART 2: Using a moving average on the rectified, reconstructed EMG signal to get an envelope
Rect_recoEMG = abs(recoEMG); % Rectify the EMG
MA_coef_num = 5000; % Length of the moving average filter (NOTE: Length determines cut-off frequency)
MA = ones(1,MA_coef_num)/MA_coef_num; % Impulse response of the moving average filter (see slide 30)
% Plot the frequency response of the moving average filter (see slide 31).
% NOTE: The absolute value of the frequency response is plotted in log
% scale
figure(7)
freqz(MA);
% Compute the envelope of the EMG by filtering the rectified signal
[Here calculate the envelope of the EMG signal]
% Plot the envelope of the signal together with torque
[Here plot the EMG envelope with superimposed torque]


         
加QQ：99515681  WX：codinghelp  Email: 99515681@qq.com
