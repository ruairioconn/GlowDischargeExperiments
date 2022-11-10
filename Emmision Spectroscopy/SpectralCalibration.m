%% Calibrate and Save Spectral Measurements
tic
clear;clc;
calLamp = 0;
save_file = 0;

% Uncalibrated
BB_dir = "2022-09-15 - Spectra/Blackbody/"; % Blackbody directory
BB_background_dir = "2022-09-15 - Spectra/Blackbody Background/"; % Blackbody background directory
GD_dir = "2022-09-15 - Spectra/Glow Discharge/100 mTorr/5W/"; % Data directory
GD_background_dir = "2022-09-15 - Spectra/Glow Discharge Background/Light Off"; % Data background directory
BB_int_time = 1000E-3; % Integration time (s)
BB_background_int_time = 1000E-3; % Integration time (s)
GD_int_time = 15000E-3; % Integration time (s)
GD_background_int_time = 5000E-3; % Integration time (s)
BB_temp = 1200; % Blackbody temperature (C)
L = 10E-2; % Depth of plasma slab (or Depth of Field if smaller)
save_file_name = "2022-09-15 - Spectra/Calibrated Spectra/100 mTorr/100mTorr-5W.txt";

% Wavelength cal (from cal lamp)
p = [1.0015, -5.2910]; % From 765nm central wavelength at 300 g/mm
p = [1.0015, 0.2164]; % From Step and Glue 650-850nm at 300 g/mm
p = [0.9993, 1.7884]; % From Step and Glue 650-950nm at 300 g/mm
p = [0.9966 4.2216]; % For 2022-09-15 Spectra
if calLamp == 1
    Lamp_dir = "2022-09-15 - Spectra/Cal Lamp/";
    comp_Lamp = spec_image_comp(Lamp_dir); % Compile images from directory
    [Lamp_lambda, Lamp_intensity] = spec_array2list(comp_Lamp);
    [NIST_lambda_cal, Lamp_lambda_cal, NIST_ix_cal, data_ix_cal, diff_cal] = chooseTransitions(Lamp_lambda,Lamp_intensity,"rough");
    close all
    p = polyfit(Lamp_lambda_cal+diff_cal,NIST_lambda_cal,1);
end

%% Load in data

comp_BB = spec_image_comp(BB_dir); % Compile images from directory
comp_BB(:,4) = comp_BB(:,4)/BB_int_time; % Normalize by integration time
comp_BB_background = spec_image_comp(BB_background_dir); % Compile images from directory
comp_BB_background(:,4) = comp_BB_background(:,4)/BB_background_int_time; % Normalize by integration time
BB_correct = spec_array_subtract(comp_BB,comp_BB_background); % Subtract background from data
[BB_lambda, BB_intensity] = spec_array2list(BB_correct); % Convert image array into lists


comp_GD = spec_image_comp(GD_dir); % Compile images from directory
comp_GD(:,4) = comp_GD(:,4)/GD_int_time; % Normalize by integration time
comp_GD_background = spec_image_comp(GD_background_dir); % Compile images from directory
comp_GD_background(:,4) = comp_GD_background(:,4)/GD_background_int_time; % Normalize by integration time
GD_correct = spec_array_subtract(comp_GD,comp_GD_background); % Subtract background from data
[GD_lambda, GD_intensity] = spec_array2list(GD_correct); % Convert image aray into lists
[GD_intensity,B,baseline] = BaselineCorrect(GD_lambda,GD_intensity);
BB_intensity = BB_intensity - baseline;
% GD_intensity = detrend(detrend(GD_intensity,2),2); % Remove baseline from data
% GD_intensity = GD_intensity - min(GD_intensity);
% GD_intensity = detrend(GD_intensity);

% GD_lambda: uncalibrated wavelength of data
% GD_intensity: uncalibrated intensity of data
% BB_lambda: uncalibrated wavelength of blackbody data
% BB_intensity: uncalibrated intensity of blackbody data

%% Calibration
% Apply linear fit to wavelength
lambda = polyval(p,GD_lambda);

% Apply blackbody correction to data
[calibrated_BB_intensity, BB_theoretical_intensity, scale_factor] = spec_blackbody_correct(lambda,BB_intensity,BB_intensity, BB_temp);
[calibrated_intensity, BB_theoretical_intensity, data_scale_factor] = spec_blackbody_correct(lambda,GD_intensity,BB_intensity, BB_temp);

% lambda: calibrated wavelength of data
% calibrated_intensity: calibrated intensity of data

%% Save to file

if save_file == 1
    T = table(lambda',calibrated_intensity');
    T.Properties.VariableNames = {'Wavelength','Intensity'};
    writetable(T,save_file_name);
end

toc