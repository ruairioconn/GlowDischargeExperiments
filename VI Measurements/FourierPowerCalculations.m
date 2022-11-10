%% MATLAB I-V Data Processing

clear;clc;

% Experiment details
channels = {'V_m','I_m','dV','dI'};
powers = ["5W";"10W";"15W";"20W";"25W";"30W";"35W";"40W";"45W";"50"];
powers_in = [5,10,15,20,25,30,35,40,45,50];
nRuns = 4;
fs = 1E9;
nSamples = 10000;
N_harm = 5;
pressure = "250mTorr";
frequency = 13.56E6;
phase_err = 0.7265;
phase_err = 57.887*pi/180;

% Folder path to data
files = "C:\Users\Ruairi O'Connor\Box\Glow Discharge Data\IV Data\2022-09-15\" + pressure + "\*.txt";
fileList = dir(files);
number_runs = length(fileList)/length(powers);

%Arrays to populate
V_RMS_raw = [];
V_P2P_raw = [];
V_DC_raw = [];
freq_raw = [];
phase_raw = [];
Power_raw = [];

FourierAmp_raw = [];
FourierPhase_raw = [];
FourierFreq_raw = [];

%Read data for each run into arrays
for i=1:1:length(fileList)
    fileName = append(fileList(i).folder,'\',fileList(i).name);
    data = readtable(fileName);

    V = data.(1);
    I = data.(2);
    dV = data.(3);
    dI = data.(4);
    
    [freqI,ampI,phaseI] = harmonicAnalysis(I,frequency,N_harm,fs);
    [freqV,ampV,phaseV] = harmonicAnalysis(V,frequency,N_harm,fs);
    
    phase_calc_raw(i,:) = phaseI - phaseV + phase_err;

    V_DC_raw(i) = mean(V);
    V_P2P_raw(i) = peak2peak(V);
    I_P2P_raw(i) = peak2peak(I);
    V_RMS_raw(i) = rms(V);
    I_RMS_raw(i) = rms(I);
    P_Fourier_raw(i) = sum(sqrt(2).*ampV.*ampI.*cos(phase_calc_raw(i)));

end

%Read averaged at each power into arrays
for i = 1:1:length(powers)
    % (i-1)*number_runs+1:i*number_runs gives 1,2,3,4;5,6,7,8;9,10,11,12;etc..
    V_RMS(i,:) = mean(V_RMS_raw((i-1)*number_runs+1:i*number_runs));
    V_P2P(i,:) = mean(V_P2P_raw((i-1)*number_runs+1:i*number_runs));
    I_RMS(i,:) = mean(I_RMS_raw((i-1)*number_runs+1:i*number_runs));
    I_P2P(i,:) = mean(I_P2P_raw((i-1)*number_runs+1:i*number_runs));
    V_DC(i,:) = mean(V_DC_raw((i-1)*number_runs+1:i*number_runs));
    phase_calc(i,:) = mean(phase_calc_raw((i-1)*number_runs+1:i*number_runs));
    P_Fourier(i,:) = mean(P_Fourier_raw((i-1)*number_runs+1:i*number_runs));
end

% figure
% scatter(V_P2P,abs(P_Fourier))
% xlabel('Electrode Voltage - P2P (V)')
% ylabel('Measured Power (W)')
% title('Plamsa Power: ' + pressure + ', 13.56 MHz (uncorrected)')

%% CLC Correction

%System parameters from VNA fit
C_sys = 114E-12;
L_loop = 250E-9;
C_coax = 48.9E-12;
R_loop = 2.32;
R_coax = 1.46;
L_coax = 147E-9;
x = [C_sys,L_loop,C_coax,R_loop,R_coax,L_coax];

% Component value uncertainties
dC_sys = 1.212E-12;
dL_loop = 1.852E-9;
dC_coax = 3.431E-13;
dR_loop = 0.1326;
dR_coax = 0.2387;
dL_coax = 2.4212E-9;
std = [dC_sys,dL_loop,dC_coax,dR_loop,dR_coax,dL_coax];

% Calculate corrected values
Vm = 0.5.*V_P2P.*exp(1i.*0);
Im = 0.5.*I_P2P.*exp(1i.*phase_calc);
[Vp,Ip,Vp_unc,Ip_unc] = VI(Vm,Im,x,std,13.56E6);
Pp = 0.5.*abs(Vp).*abs(Ip).*cos(-angle(Vp)+angle(Ip));

dPpdVp = (cos(angle(Ip) - angle(Vp)).*abs(Ip).*sign(Vp))/2 - (sin(angle(Ip) - angle(Vp)).*imag(Vp).*abs(Ip).*abs(Vp))./(2.*(real(Vp).^2 + imag(Vp).^2));
dPpdIp = (cos(angle(Ip) - angle(Vp)).*abs(Vp).*sign(Ip))/2 + (imag(Ip).*sin(angle(Ip) - angle(Vp)).*abs(Ip).*abs(Vp))./(2.*(imag(Ip).^2 + real(Ip).^2));
Pp_unc = sqrt((dPpdVp.*Vp_unc).^2+(dPpdIp.*Ip_unc).^2);

yneg = abs(Vp_unc);
ypos = abs(Vp_unc);
xneg = abs(Ip_unc);
xpos = abs(Ip_unc);

% Plot results
figure
errorbar(abs(Ip),abs(Vp),yneg,ypos,xneg,xpos,'-o')
title('I-V (1 Torr)')
xlabel('Current (A)')
ylabel('Voltage (V)')

yneg = abs(Pp_unc);
ypos = abs(Pp_unc);
xneg = abs(Vp_unc);
xpos = abs(Vp_unc);

figure
errorbar(abs(Vp),abs(Pp),yneg,ypos,xneg,xpos,'-o')
title('Plasma Power (1 Torr)')
xlabel('Voltage (V)')
ylabel('Power (W)')

P = polyfit(abs(Ip),abs(Vp),1);
R = P(1)

function [Vp,Ip, Vp_unc,Ip_unc] = VI(Vm,Im,x,std,freq)
w = 2.*pi.*freq;
C_sys = x(1);
L_loop = x(2);
C_coax = x(3);
R_loop = x(4);
R_coax = x(5);
L_coax = x(6);

% Circuit model correction
I_coax = Im - (Vm-1i.*w.*L_coax.*Im)./(R_coax - (1./(1i.*w.*C_coax)));
Vp = Vm - I_coax.*(R_loop+1i.*w.*L_loop)-1i.*w.*L_coax.*Im;
Ip = Im - I_coax + 1i.*w.*C_sys.*(-Vm+I_coax.*(R_loop+1i.*w.*L_loop)+1i.*w.*L_coax.*Im);

dC_sys = std(1);
dL_loop = std(2);
dC_coax = std(3);
dR_loop = std(4);
dR_coax = std(5);
dL_coax = std(6);
dVm = abs(Vm)*0.03;
dIm = abs(Im)*0.01;

% Uncertainty characterization

dVpdC_sys = 0;
dVpdL_loop = -w*(Im - (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w)))*1i;
dVpdC_coax = ((Vm - Im*L_coax*w*1i)*(R_loop + L_loop*w*1i)*1i)/(C_coax^2*w*(R_coax + 1i/(C_coax*w))^2);
dVpdR_loop = - Im + (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w));
dVpdR_coax = -((Vm - Im*L_coax*w*1i)*(R_loop + L_loop*w*1i))/(R_coax + 1i/(C_coax*w))^2;
dVpdL_coax = - Im*w*1i - (Im*w*(R_loop + L_loop*w*1i)*1i)/(R_coax + 1i/(C_coax*w));
dVpdVm = (R_loop + L_loop*w*1i)/(R_coax + 1i/(C_coax*w)) + 1;
dVpdIm = - L_coax*w*1i - ((L_coax*w*1i)/(R_coax + 1i/(C_coax*w)) + 1)*(R_loop + L_loop*w*1i);

Vp_unc = sqrt((dVpdC_sys*dC_sys).^2+(dVpdL_loop*dL_loop).^2+(dVpdC_coax*dC_coax).^2+(dVpdR_loop*dR_loop).^2+(dVpdR_coax*dR_coax).^2+(dVpdL_coax*dL_coax).^2+(dVpdVm*dVm).^2+(dVpdIm*dIm).^2);

dIpdC_sys = w*(- Vm + (Im - (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w)))*(R_loop + L_loop*w*1i) + Im*L_coax*w*1i)*1i;
dIpdL_loop = -C_sys*w^2*(Im - (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w)));
dIpdC_coax = ((Vm - Im*L_coax*w*1i)*1i)/(C_coax^2*w*(R_coax + 1i/(C_coax*w))^2) + (C_sys*(Vm - Im*L_coax*w*1i)*(R_loop + L_loop*w*1i))/(C_coax^2*(R_coax + 1i/(C_coax*w))^2);
dIpdR_loop = C_sys*w*(Im - (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w)))*1i;
dIpdR_coax = - (Vm - Im*L_coax*w*1i)/(R_coax + 1i/(C_coax*w))^2 + (C_sys*w*(Vm - Im*L_coax*w*1i)*(R_loop + L_loop*w*1i)*1i)/(R_coax + 1i/(C_coax*w))^2;
dIpdL_coax = - (Im*w*1i)/(R_coax + 1i/(C_coax*w)) + C_sys*w*(Im*w*1i + (Im*w*(R_loop + L_loop*w*1i)*1i)/(R_coax + 1i/(C_coax*w)))*1i;
dIpdVm = 1/(R_coax + 1i/(C_coax*w)) - C_sys*w*((R_loop + L_loop*w*1i)/(R_coax + 1i/(C_coax*w)) + 1)*1i;
dIpdIm = - (L_coax*w*1i)/(R_coax + 1i/(C_coax*w)) + C_sys*w*(L_coax*w*1i + ((L_coax*w*1i)/(R_coax + 1i/(C_coax*w)) + 1)*(R_loop + L_loop*w*1i))*1i;

Ip_unc = sqrt((dIpdC_sys*dC_sys).^2+(dIpdL_loop*dL_loop).^2+(dIpdC_coax*dC_coax).^2+(dIpdR_loop*dR_loop).^2+(dIpdR_coax*dR_coax).^2+(dIpdL_coax*dL_coax).^2+(dIpdVm*dVm).^2+(dIpdIm*dIm).^2);

end
