%% Calculate Parameters from Emission Spectra
clear;clc;

choose_transitions = 0;
data = readtable("2022-09-15 - Spectra/Calibrated Spectra/1 Torr/1Torr-50W.txt");
lambda = data.Wavelength;
calibrated_intensity = data.Intensity;
L = 32.5E-3;
NIST_data = call_NIST(lambda,"Table");

if choose_transitions == 0
    NIST_lambda = [727.2936, 731.1716, 735.3293, 737.2118, 738.3980, 739.2980, 743.6297, 750.3869, 751.4652, 763.5106, 772.3761, 794.8176, 800.6157, 801.4786];
    NIST_lambda = [667.7282  675.6163  687.1289  696.5431  703.0251  706.7218  714.7042  720.6980  726.5172  727.2936  731.1716  735.3293  737.2118  738.3980  743.5368  750.3869  751.4652  763.5106 772.3761  789.1075  794.8176  800.6157  801.4786]';
    Ar_2s4_transitions = [696.5431,706.7218,714.7042,763.5106,772.3761,801.4786,811.5311,912.2967];
    Ar_2s3_transitions = [667.7282,727.2936,738.3980,751.4652,800.6157,810.3693,842.4648,965.7784];
    Ar_2s2_transitions = [772.4207,794.8176,866.7944];
    Ar_2s1_transitions = [750.3869,826.4522,840.8210,852.1442,922.4499,935.4220,978.5403];
    Ar_4p_4s_transitions = sort([Ar_2s1_transitions,Ar_2s2_transitions,Ar_2s3_transitions,Ar_2s4_transitions]);
    NIST_lambda = Ar_4p_4s_transitions(Ar_4p_4s_transitions>=min(lambda) & Ar_4p_4s_transitions<=max(lambda));
else
    [NIST_lambda, data_lambda, NIST_ix, data_ix, diff, pks, locs, wdths, A] = chooseTransitions(lambda',calibrated_intensity');
end

%% Spectral Fitting

% [pks,locs,wdths] = findpeaks(calibrated_intensity,lambda,'MinPeakHeight',0.01*max(calibrated_intensity),'MinPeakWidth',0.2);
[pks,locs,wdths] = findpeaks(calibrated_intensity,lambda,'MinPeakWidth',0.3);

NIST_ix = [];
data_ix = [];
for i = 1:length(NIST_lambda)
    NIST_ix = [NIST_ix, find(NIST_data.lambda==NIST_lambda(i))];
    [~,ixmin] = min(abs(locs-NIST_lambda(i)));
    data_ix = [data_ix, ixmin];
end
A = NIST_data.A_ki(NIST_ix)';
data_lambda = locs(data_ix);
pks = pks(data_ix);
locs = locs(data_ix);
wdths = wdths(data_ix);

% Select transitions for spectral fitting and collect initial solver guesses
% [NIST_lambda, data_lambda, NIST_ix, data_ix, diff, pks, locs, wdths, A] = chooseTransitions(lambda,calibrated_intensity);
% close all
initguess = [pks,locs,wdths/2,wdths/2];
initguess = [pks,locs,wdths];
initguess(22,2) = 911;

% Fit Gaussian to spectrum from initial guesses
[spectrum, fit_params, res, I] = SpectralFitting(lambda, calibrated_intensity, initguess, "Gauss");
n = length(pks); % number of transitions
a = fit_params(1:n); % transition peak height
p0 = fit_params(n+1:2*n); % center wavelength
w = fit_params(2*n+1:3*n); % width of line

% % Calculate line intensities
% for i = 1:length(pks)
%     I(i) = (sqrt(pi)*a(i)*(w(i)*10^-9))/(2*sqrt(log(2)));
% end

% I: intensity of each transition line (area under curve)
% spectrum: fitted Gaussian result
% NIST_lambda: wavelengths of chosen transitions
% NIST_ix: indices of chosen transitions in NIST table
% A: Einstein A coefficients for chosen transitions

figure
hold on
plot(lambda,calibrated_intensity)
plot(lambda,spectrum)
xlabel('Wavelength (nm)')
ylabel('Intensity (W/m3-sr)')
legend('Measurement (calibrated)','Fit')

%% Calculations

% Calculate number densities
NIST_out = call_NIST(lambda,NIST_lambda);
h = 6.626E-34;
c = 299.8E6;
n_dens = (4*pi*I.*(p0*10^-9))./(A.*h.*c.*L);
n_dens_sum = [mean(n_dens([22])),mean(n_dens([16])),mean(n_dens([14,19])),mean(n_dens([10,15,21,24])),mean(n_dens([9,13,23])),mean(n_dens([8])),mean(n_dens([4,12,20])),mean(n_dens([3,6,18])),mean(n_dens([2,5,11,17])),mean(n_dens([1,7]))];
n_dens_sum = [mean(n_dens([22])),mean(n_dens([16])),mean(n_dens([14,19])),mean(n_dens([10,15,21])),mean(n_dens([9,13])),mean(n_dens([8])),mean(n_dens([12,20])),mean(n_dens([3,6,18])),mean(n_dens([2,5,11,17])),mean(n_dens([7]))];

n_dens_output = [n_dens_sum,sum(n_dens_sum)]
sum(n_dens_sum)

% Calculate uncertainty
dndA = -(4*pi*I.*(p0*10^-9))./(A.^2.*h.*c.*L);
dndI = (4*pi.*(p0*10^-9))./(A.*h.*c.*L);
dndlambda = (4*pi*I)./(A.*h.*c.*L);
dA = NIST_out.accuracy.*A;
dI = zeros(size(dA));
dlambda = zeros(size(dA));

n_dens_uncertainty = sqrt((dndA.*dA).^2+(dndI.*dI).^2+(dndlambda.*dlambda).^2);

% Boltzmann plot
g = NIST_out.g;
E = NIST_out.E_upper;
Boltz_y = log(I.*(NIST_lambda*10^-9)./(A.*g'));
% 
% good = [1,2,3,6,8,9,10];
% bad = [4,5,7];

p = polyfit(E,Boltz_y,1);
L = polyval(p,E);

k = 8.617E-5;
T = -1/(k*p(1));

% figure
% hold on
% scatter(E,Boltz_y,'blue','filled')
% % scatter(E(bad),Boltz_y(bad),'red','filled')
% % labelpoints(E,Boltz_y,NIST_lambda,'N',0.1)
% plot(E,L)

% NIST data for wavelength range
NIST_data = call_NIST(lambda,"data");

% n_dens: number density corresponding to each transition chosen
% NIST_data: table of NIST data for all transitions in wavelength range
%       e.g.: NIST_data.lambda(NIST_ix): wavelengths of all transitions
%       chosen

%% Metastable Densities

% lambda = NIST_data.lambda([126 161]);
% A = NIST_data.A_ki([126 161]);
% g_i = NIST_data.g_i([126 161]);
% g_k = NIST_data.g_k([126 161]);
% input = [lambda(:),A(:),g_k(:),g_i(:)];
% 
% n0 = [1E15,1E15];
% ratio0 = pks(18)/pks(22);
% 
% fun = @metastableDens;
% ratio = fun(n0,input);
% 
% opts = statset('nlinfit');
% opts.Method = 'NonlinearLeastSquares';
% beta = nlinfit(input,ratio0,fun,n0,opts);
% 
% function ratio = metastableDens(n, input)
% 
% lambda = input(:,1);
% A_ij = input(:,2);
% g_i = input(:,3);
% g_j = input(:,4);
% 
% T = 300;
% m = 6.6335E-26;
% k = 1.38E-26;
% l = 10E-2;
% 
% % Transition 1
% P_ij = lambda(1)*sqrt(m/(2*pi*k*T));
% K_ij = ((lambda(1)^2)/(8*pi))*P_ij*(g_i(1)/g_j(1))*n(1)*A_ij(1);
% gamma_ij = (2-exp((-10^-3)*K_ij*l))/(1+K_ij*lambda(1));
% Phi_ij = gamma_ij*A_ij(1);
% 
% % Transition 1
% P_ik = lambda(2)*sqrt(m/(2*pi*k*T));
% K_ik = ((lambda(2)^2)/(8*pi))*P_ik*(g_i(2)/g_j(2))*n(2)*A_ij(2);
% gamma_ik = (2-exp((-10^-3)*K_ik*l))/(1+K_ik*lambda(2));
% Phi_ik = gamma_ik*A_ij(2);
% 
% ratio = Phi_ij/Phi_ik;
% 
% end
