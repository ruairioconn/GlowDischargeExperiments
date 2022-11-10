%% Absorption Estimation

clear;clc;

data = readtable("2022-09-15 - Spectra/Calibrated Spectra/1 Torr/1Torr-5W.txt");
lambda = data.Wavelength;
% lambda = linspace(min(lambda_data),max(lambda_data),100000);
calibrated_intensity = data.Intensity;
NIST_data = call_NIST(lambda,"Table");

Ar_2s4_transitions = [696.5431,706.7218,714.7042,763.5106,772.3761,801.4786,811.5311,912.2967];
Ar_2s3_transitions = [667.7282,727.2936,738.3980,751.4652,800.6157,810.3693,842.4648,965.7784];
Ar_2s2_transitions = [772.4207,794.8176,866.7944];
Ar_2s1_transitions = [750.3869,826.4522,840.8210,852.1442,922.4499,935.4220,978.5403];
Ar_4p_4s_transitions = sort([Ar_2s1_transitions,Ar_2s2_transitions,Ar_2s3_transitions,Ar_2s4_transitions]);
NIST_lambda = Ar_4p_4s_transitions(Ar_4p_4s_transitions>=min(lambda) & Ar_4p_4s_transitions<=max(lambda));

states = {"2p10";"2p9";"2p8";"2p7";"2p6";"2p5";"2p4";"2p3";"2p2";"2p1"};
n_dens_4p = [6.485E12,6.249E12,6.8774E12,6.819E12,8.884E12,2.014E12,3.00233E12,4.61467E12,3.276E12,3.243E12];
n_dens_4s = [3.26E16,6.20E15,6.52E15,6.20E15];

Ar2p10 = [912.2967;n_dens_4p(1);n_dens_4s(1)];
Ar2p9 = [811.5311;n_dens_4p(2);n_dens_4s(1)];
Ar2p8 = [801.4786,842.4648;n_dens_4p(3),n_dens_4p(3);n_dens_4s(1),n_dens_4s(2)];
Ar2p7 = [772.3761,810.3693,866.7944,935.4220;n_dens_4p(4),n_dens_4p(4),n_dens_4p(4),n_dens_4p(4);n_dens_4s(1),n_dens_4s(2),n_dens_4s(3),n_dens_4s(4)];
Ar2p6 = [763.5106,800.6157,922.4499;n_dens_4p(5),n_dens_4p(5),n_dens_4p(5);n_dens_4s(1),n_dens_4s(2),n_dens_4s(4)];
Ar2p5 = [751.4652;n_dens_4p(6);n_dens_4s(2)];
Ar2p4 = [714.7042,794.8176,852.1442;n_dens_4p(7),n_dens_4p(7),n_dens_4p(7);n_dens_4s(1),n_dens_4s(3),n_dens_4s(4)];
Ar2p3 = [706.7218,738.3980,840.8210;n_dens_4p(8),n_dens_4p(8),n_dens_4p(8);n_dens_4s(1),n_dens_4s(2),n_dens_4s(4)];
Ar2p2 = [696.5431,727.2936,772.4207,826.4522;n_dens_4p(9),n_dens_4p(9),n_dens_4p(9),n_dens_4p(9);n_dens_4s(1),n_dens_4s(2),n_dens_4s(3),n_dens_4s(4)];
Ar2p1 = [667.7282,750.3869;n_dens_4p(10),n_dens_4p(10);n_dens_4s(2),n_dens_4s(4)];

% Ar4p = table(num2str(Ar2p10),num2str(Ar2p9),num2str(Ar2p8),num2str(Ar2p7),num2str(Ar2p6),num2str(Ar2p5),num2str(Ar2p4),num2str(Ar2p3),num2str(Ar2p2),num2str(Ar2p1));
Ar4p = table(Ar2p10,Ar2p9,Ar2p8,Ar2p7,Ar2p6,Ar2p5,Ar2p4,Ar2p3,Ar2p2,Ar2p1);

m = 6.634E-26; % Ar Mass (kg)
% m = 39.948; % Ar Mass (amu)
h = 4.135667696E-15; % Planck's Constant (eV/Hz)
% h = 6.626E-34; % Planck's Constant (m2.kg/s)
k = 1.38E-23; % Boltzmann's Constant (m2.kg/s2.K)
c = 299792458; % Speed of Light (m/s)
T = 10000*300; % Bulk Gas Temperature (K)
L = 5E-2; % Path Length (m)

I_abs = zeros(size(lambda));
nu = c./(lambda*1E-9);

% for i=1:length(NIST_lambda)
%     NIST_ix = find(NIST_data.lambda==NIST_lambda(i));
%     E_l = NIST_data.E_i(NIST_ix);
%     E_u = NIST_data.E_k(NIST_ix);
%     g_l = NIST_data.g_i(NIST_ix);
%     g_u = NIST_data.g_k(NIST_ix);
%     A_ul = NIST_data.A_ki(NIST_ix);
% 
%     index = find(Ar4p{1,:}==NIST_lambda(i));
%     n_l = Ar4p{3,:}(index);
%     n_u = Ar4p{2,:}(index);
%     n_abs = n_l;
% 
% %     nu0 = (E_u-E_l)/h;
% %     delta_nu_D = nu0*sqrt(2*k*T/m)/c;
% % %     delta_nu_D = nu0*(4.103E-7)*sqrt(T/m);
% %     Phi = (1/(sqrt(pi)*delta_nu_D))*exp(-((nu-nu0)./(delta_nu_D)).^2);
%     lambda0 = NIST_lambda(i)*1E-9;
%     B_lu = (g_u/g_l)*((c^3)/(8*pi*h*(c/lambda0)^3))*A_ul;
% %     Phi = (c/lambda0)*sqrt(m/(2*pi*k*T))*exp(-m*(c*(1-lambda*1E-9./lambda0)).^2./(2*k*T));
% %     S_lu = (n_l/n_abs)*(1-((n_u*g_l)/(n_l*g_u)))*(g_u/g_l)*(c^2/(8*pi*nu0^2))*A_ul;
%     delta_lambda_D = lambda0*2*sqrt(2*log(2)*((k*T)/(m*c^2)));
%     delta_lambda_D = 0.37E-9;
%     Phi = (1/(sqrt(pi)*delta_lambda_D))*exp(-((lambda*1E-9 - lambda0)./delta_lambda_D).^2);
% %     S_lu = (n_l/n_abs)*h*(c/lambda0)*(1-((n_u*g_l)/(n_l*g_u)))*B_lu;
% %     alpha = n_abs*S_lu*Phi;
%     alpha = (h/lambda0)*B_lu*(n_l-n_u*(g_l/g_u))*Phi;
% 
%     Tau = exp(-alpha*L);
%     I_abs = I_abs + (Tau-1);
% end

Phi = zeros(size(lambda));

for i=1:length(NIST_lambda)
    NIST_ix = find(NIST_data.lambda==NIST_lambda(i));
    E_l = NIST_data.E_i(NIST_ix);
    E_u = NIST_data.E_k(NIST_ix);
    g_l = NIST_data.g_i(NIST_ix);
    g_u = NIST_data.g_k(NIST_ix);
    A_ul = NIST_data.A_ki(NIST_ix);

    index = find(Ar4p{1,:}==NIST_lambda(i));
    n_l = Ar4p{3,:}(index);
    n_u = Ar4p{2,:}(index);
    n_abs = n_l;

    nu0 = (E_u-E_l)/h;
    delta_nu_D = nu0*sqrt(2*k*T/m)/c;
%     delta_nu_D = nu0*(4.103E-7)*sqrt(T/m);
    Phi = Phi + (1/(sqrt(pi)*delta_nu_D))*exp(-((nu-nu0)./(delta_nu_D)).^2);
    S_lu = (n_l/n_abs)*(1-((n_u*g_l)/(n_l*g_u)))*(g_u/g_l)*(c^2/(8*pi*nu0^2))*A_ul;
    alpha = n_abs*S_lu*Phi;

    Tau = exp(-alpha*L);
    I_abs = I_abs + (Tau-1);
end

I_abs = I_abs + 1;

figure
plot(lambda,I_abs)

figure
plot(lambda,calibrated_intensity)
hold on
plot(lambda,calibrated_intensity.*I_abs)

p = 1 * 133.3;
n0 = p/(k*T);
k_eV = 8.617E-5;
h = 6.626E-34;
T_el = 0.6*11604;
wdth = 0.3352;
L = 10E-2;

S_model = zeros(size(lambda));

for i = 1:length(NIST_lambda)
    NIST_ix = find(NIST_data.lambda==NIST_lambda(i));
    E_u = NIST_data.E_k(NIST_ix);
    A = NIST_data.A_ki(NIST_ix);
    n_4p(i) = n0*exp(-E_u/(k_eV*T_el));

    G0 = [1,NIST_lambda(i),wdth];
    S_n_profile = (2.*sqrt(log(2))./(sqrt(pi).*wdth.*1E-9)).*(n_4p(i)).*(A./(4.*pi)).*((h.*c)./(NIST_lambda(i).*1E-9)).*Gauss(G0,lambda);
    S_n = S_n_profile*L;
    S_model = S_model + S_n;
end

figure
plot(lambda,calibrated_intensity/calibrated_intensity(1356))
hold on
plot(lambda,S_model/S_model(1355))
legend("Measurement","Theory")
