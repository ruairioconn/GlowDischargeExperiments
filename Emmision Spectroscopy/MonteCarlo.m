%% Monte Carlo Uncertainty of Ar (4p) Population Densities

clear;clc;
tic

density_input = readmatrix('2022-09-15 - Spectra\250 mTorr densities.csv')*1E13;
write_file = '2022-09-15 - Spectra\MC Results\250 mTorr\5W_uncertainties.txt';
n_dens = [6.485E12,6.249E12,6.8774E12,6.819E12,8.884E12,2.014E12,3.00233E12,4.61467E12,3.276E12,3.243E12];
n_dens = [0.0384, 2.2100, 2.1030 ,1.2869, 1.6022, 0.6411, 2.1157, 0.7266, 0.7611, 3.6023]*1E12;
n_dens = density_input(1,1:10);

numTrials = 10000;

id = 'instrument:udp:ClassToBeRemoved';
warning('off',id)
id = 'MATLAB:table:ModifiedAndSavedVarnames';
warning('off',id)

lambda = linspace(650,950,1024*3);

include = ["A","L","Cal"];

L = 32.5E-3;
NIST_data = call_NIST(lambda,"Table");
h = 6.626E-34;
c = 299.8E6;
k_B = 1.38E-23;
wdth = 0.3352;
t = 5000E-3;
A_px = (12.8E-6)^2;
eta_csv = readtable("PiMax4 QE.csv");
dL = 8.6E-3;

eta_lambda = [eta_csv.(1); linspace(900.5,1000,100)'];
% eta = [eta_csv.(2); (eta_csv.(2)(end)*exp(eta_csv.(1)(end)/10))*exp(-linspace(900.5,1000,100)'/10)]/100;
eta = [eta_csv.(2); ones(100,1)]/100;

Ar_2s4_transitions = [696.5431,706.7218,714.7042,763.5106,772.3761,801.4786,811.5311,912.2967];
Ar_2s3_transitions = [667.7282,727.2936,738.3980,751.4652,800.6157,810.3693,842.4648,965.7784];
Ar_2s2_transitions = [772.4207,794.8176,866.7944];
Ar_2s1_transitions = [750.3869,826.4522,840.8210,852.1442,922.4499,935.4220,978.5403];
Ar_4p_4s_transitions = sort([Ar_2s1_transitions,Ar_2s2_transitions,Ar_2s3_transitions,Ar_2s4_transitions]);
NIST_lambda = Ar_4p_4s_transitions(Ar_4p_4s_transitions>=min(lambda) & Ar_4p_4s_transitions<=max(lambda));

states = {"2p_1_0";"2p_9";"2p_8";"2p_7";"2p_6";"2p_5";"2p_4";"2p_3";"2p_2";"2p_1"};

Ar2p10 = [912.2967;n_dens(1)];
Ar2p9 = [811.5311;n_dens(2)];
Ar2p8 = [801.4786,842.4648;n_dens(3),n_dens(3)];
Ar2p7 = [772.3761,810.3693,866.7944,935.4220;n_dens(4),n_dens(4),n_dens(4),n_dens(4)];
Ar2p6 = [763.5106,800.6157,922.4499;n_dens(5),n_dens(5),n_dens(5)];
Ar2p5 = [751.4652;n_dens(6)];
Ar2p4 = [714.7042,794.8176,852.1442;n_dens(7),n_dens(7),n_dens(7)];
Ar2p3 = [706.7218,738.3980,840.8210;n_dens(8),n_dens(8),n_dens(8)];
Ar2p2 = [696.5431,727.2936,772.4207,826.4522;n_dens(9),n_dens(9),n_dens(9),n_dens(9)];
Ar2p1 = [667.7282,750.3869;n_dens(10),n_dens(10)];

Ar4p = table(Ar2p10,Ar2p9,Ar2p8,Ar2p7,Ar2p6,Ar2p5,Ar2p4,Ar2p3,Ar2p2,Ar2p1);

NIST_ix = [];
for i = 1:length(NIST_lambda)
    NIST_ix = [NIST_ix, find(NIST_data.lambda==NIST_lambda(i))];
end
A = NIST_data.A_ki(NIST_ix)';
NIST_out = call_NIST(lambda,NIST_lambda);
A_accuracy = NIST_out.accuracy;
g = NIST_data.g_k(NIST_ix);
E = NIST_data.E_k(NIST_ix);

n_dens_MC = zeros(numTrials,25);
profile = "Uniform";
true_shape = "Gauss";
assumed_shape = "Gauss";

display("Main Monte Carlo Analysis")

ppm = ParforProgressbar(numTrials);
parfor i=1:numTrials
    n_dens_MC(i,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile,t,A_px,eta_lambda,eta,dL,include,true_shape,assumed_shape);
    ppm.increment();
end
delete(ppm);

MC_mean = mean(n_dens_MC,1);
MC_std = std(n_dens_MC,0,1);
scatter_x = [1,2,3,3,4,4,4,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,10,10];
scatter_y = MC_mean([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7]);
n_dens_mean = [mean(MC_mean([22])),mean(MC_mean([16])),mean(MC_mean([14,19])),mean(MC_mean([10,15,21,24])),mean(MC_mean([9,13,23])),mean(MC_mean([8])),mean(MC_mean([4,12,20])),mean(MC_mean([3,6,18])),mean(MC_mean([2,5,11,17])),mean(MC_mean([1,7]))];
n_dens_std = [mean(MC_std([22])),mean(MC_std([16])),mean(MC_std([14,19])),mean(MC_std([10,15,21,24])),mean(MC_std([9,13,23])),mean(MC_std([8])),mean(MC_std([4,12,20])),mean(MC_std([3,6,18])),mean(MC_std([2,5,11,17])),mean(MC_std([1,7]))];

writematrix(n_dens_MC, write_file);



%% Lineshape Effects

% line_shape = ["Gauss","Lorentz","Voigt"];
%
% display("Lineshapes")
%
% ppm = ParforProgressbar(9*numTrials);
% for i = 1:3
%     display("Real shape: ", line_shape(i))
%     for j = 1:3
%         parfor k=1:numTrials
%             n_dens_MC_shape(k,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile,t,A_px,eta_lambda,eta,DoF,include,line_shape(i),line_shape(j));
%             ppm.increment();
%         end
%         MC_mean_shape = mean(n_dens_MC_shape,1);
%         MC_std_shape = std(n_dens_MC_shape,0,1);
%         n_dens_mean_shape = [mean(MC_mean_shape([22])),mean(MC_mean_shape([16])),mean(MC_mean_shape([14,19])),mean(MC_mean_shape([10,15,21,24])),mean(MC_mean_shape([9,13,23])),mean(MC_mean_shape([8])),mean(MC_mean_shape([4,12,20])),mean(MC_mean_shape([3,6,18])),mean(MC_mean_shape([2,5,11,17])),mean(MC_mean_shape([1,7]))];
%         n_dens_std_shape = [mean(MC_std_shape([22])),mean(MC_std_shape([16])),mean(MC_std_shape([14,19])),mean(MC_std_shape([10,15,21,24])),mean(MC_std_shape([9,13,23])),mean(MC_std_shape([8])),mean(MC_std_shape([4,12,20])),mean(MC_std_shape([3,6,18])),mean(MC_std_shape([2,5,11,17])),mean(MC_std_shape([1,7]))];
%         shape_n(i,j) = sum(n_dens_mean_shape);
%         shape_std(i,j) = mean(n_dens_std_shape);
%     end
% end
% delete(ppm);
%
% %% Profile Effects
%
% profiles = ["Uniform","Gaussian","Bimodal"];
%
% display("Profile Effects")
%
% ppm = ParforProgressbar(3*numTrials);
% for i = 1:3
%     display("Profile: ", profiles(i))
%         parfor k=1:numTrials
%             n_dens_MC_profile(k,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profiles(i),t,A_px,eta_lambda,eta,DoF,include,"Gauss","Gauss");
%             ppm.increment();
%         end
%         MC_mean_profile = mean(n_dens_MC_profile,1);
%         MC_std_profile = std(n_dens_MC_profile,0,1);
%         n_dens_mean_profile = [mean(MC_mean_profile([22])),mean(MC_mean_profile([16])),mean(MC_mean_profile([14,19])),mean(MC_mean_profile([10,15,21,24])),mean(MC_mean_profile([9,13,23])),mean(MC_mean_profile([8])),mean(MC_mean_profile([4,12,20])),mean(MC_mean_profile([3,6,18])),mean(MC_mean_profile([2,5,11,17])),mean(MC_mean_profile([1,7]))];
%         n_dens_std_profile = [mean(MC_std_profile([22])),mean(MC_std_profile([16])),mean(MC_std_profile([14,19])),mean(MC_std_profile([10,15,21,24])),mean(MC_std_profile([9,13,23])),mean(MC_std_profile([8])),mean(MC_std_profile([4,12,20])),mean(MC_std_profile([3,6,18])),mean(MC_std_profile([2,5,11,17])),mean(MC_std_profile([1,7]))];
%         profile_n(i) = sum(n_dens_mean_profile);
%         profile_std(i) = mean(n_dens_std_profile);
% end
% delete(ppm);


%% Error Sources
%
% numTrials_error = floor(numTrials/3);
%
% display("Error Sources")
%
% include = ["A"];
% display("A Coefficient")
% ppm = ParforProgressbar(numTrials_error);
% parfor i=1:numTrials_error
%     n_dens_MC_A(i,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile,t,A_px,eta_lambda,eta,DoF,include,true_shape,assumed_shape);
%     ppm.increment();
% end
% delete(ppm);
%
% MC_mean_A = mean(n_dens_MC_A,1);
% MC_std_A = std(n_dens_MC_A,0,1);
%
% scatter_x_A = [1,2,3,3,4,4,4,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,10,10];
% scatter_y_A = MC_mean_A([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7]);
% n_dens_mean_A = [mean(MC_mean_A([22])),mean(MC_mean_A([16])),mean(MC_mean_A([14,19])),mean(MC_mean_A([10,15,21,24])),mean(MC_mean_A([9,13,23])),mean(MC_mean_A([8])),mean(MC_mean_A([4,12,20])),mean(MC_mean_A([3,6,18])),mean(MC_mean_A([2,5,11,17])),mean(MC_mean_A([1,7]))];
% n_dens_std_A = [mean(MC_std_A([22])),mean(MC_std_A([16])),mean(MC_std_A([14,19])),mean(MC_std_A([10,15,21,24])),mean(MC_std_A([9,13,23])),mean(MC_std_A([8])),mean(MC_std_A([4,12,20])),mean(MC_std_A([3,6,18])),mean(MC_std_A([2,5,11,17])),mean(MC_std_A([1,7]))];
%
% %%%%%%%%%%%%%%%%%%%%%%%%
% % include = ["L"];
% % display("Length")
% % ppm = ParforProgressbar(numTrials_error);
% % parfor i=1:numTrials_error
% %     n_dens_MC_L(i,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile,t,A_px,eta_lambda,eta,DoF,include,true_shape,assumed_shape);
% %     ppm.increment();
% % end
% % delete(ppm);
% %
% % MC_mean_L = mean(n_dens_MC_L,1);
% % MC_std_L = std(n_dens_MC_L,0,1);
% %
% % scatter_x_L = [1,2,3,3,4,4,4,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,10,10];
% % scatter_y_L = MC_mean_L([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7]);
% % n_dens_mean_L = [mean(MC_mean_L([22])),mean(MC_mean_L([16])),mean(MC_mean_L([14,19])),mean(MC_mean_L([10,15,21,24])),mean(MC_mean_L([9,13,23])),mean(MC_mean_L([8])),mean(MC_mean_L([4,12,20])),mean(MC_mean_L([3,6,18])),mean(MC_mean_L([2,5,11,17])),mean(MC_mean_L([1,7]))];
% % n_dens_std_L = [mean(MC_std_L([22])),mean(MC_std_L([16])),mean(MC_std_L([14,19])),mean(MC_std_L([10,15,21,24])),mean(MC_std_L([9,13,23])),mean(MC_std_L([8])),mean(MC_std_L([4,12,20])),mean(MC_std_L([3,6,18])),mean(MC_std_L([2,5,11,17])),mean(MC_std_L([1,7]))];
%
% %%%%%%%%%%%%%%%%%%%%
% include = ["Cal"];
% display("Calibration Function")
% ppm = ParforProgressbar(numTrials_error);
% parfor i=1:numTrials_error
%     n_dens_MC_C(i,:) = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile,t,A_px,eta_lambda,eta,DoF,include,true_shape,assumed_shape);
%     ppm.increment();
% end
% delete(ppm);
%
% MC_mean_C = mean(n_dens_MC_C,1);
% MC_std_C = std(n_dens_MC_C,0,1);
%
% scatter_x_C = [1,2,3,3,4,4,4,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,10,10];
% scatter_y_C = MC_mean_C([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7]);
% n_dens_mean_C = [mean(MC_mean_C([22])),mean(MC_mean_C([16])),mean(MC_mean_C([14,19])),mean(MC_mean_C([10,15,21,24])),mean(MC_mean_C([9,13,23])),mean(MC_mean_C([8])),mean(MC_mean_C([4,12,20])),mean(MC_mean_C([3,6,18])),mean(MC_mean_C([2,5,11,17])),mean(MC_mean_C([1,7]))];
% n_dens_std_C = [mean(MC_std_C([22])),mean(MC_std_C([16])),mean(MC_std_C([14,19])),mean(MC_std_C([10,15,21,24])),mean(MC_std_C([9,13,23])),mean(MC_std_C([8])),mean(MC_std_C([4,12,20])),mean(MC_std_C([3,6,18])),mean(MC_std_C([2,5,11,17])),mean(MC_std_C([1,7]))];
%


% figure
% scatter([1,2,3,4,5,6,7,8,9,10],n_dens,'r*')
% hold on
% dens_scatter = scatter(scatter_x,scatter_y);
% errorbar([1,2,3,4,5,6,7,8,9,10],n_dens_mean,neg,pos,'bo','LineWidth',1)
% legend("Input","State Scatter","Monte Carlo Mean/StdDev",'FontSize',18,'Location','northwest')
% set(gca,'xtick',[1:10],'xticklabel',states,'FontSize',20)
% xlabel('Paschen State Notation','FontSize',22)
% ylabel("Population Density (#/m^3)",'FontSize',22)
% title('Monte Carlo Results','FontSize',22)
% xlim([0.9 10.1])
% % ylim([0 1.1*max(n_dens_mean+pos)])
% sadf = dataTipTextRow('Wavelength =',NIST_lambda([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7]));
% dens_scatter.DataTipTemplate.DataTipRows(end+1) = sadf;
% dens_scatter.DataTipTemplate.DataTipRows(1) = dens_scatter.DataTipTemplate.DataTipRows(1);
% dens_scatter.DataTipTemplate.DataTipRows = dens_scatter.DataTipTemplate.DataTipRows(1:2);
% dens_scatter.DataTipTemplate.DataTipRows(2).Label = 'Density = ';

% pos1_0 = 2*n_dens_std_L + 2*n_dens_std_A + 2*n_dens_std_C;
% neg1_0 = 2*n_dens_std_L + 2*n_dens_std_A + 2*n_dens_std_C;
% pos2 = 2*n_dens_std_A + 2*n_dens_std_C;
% neg2 = 2*n_dens_std_A + 2*n_dens_std_C;
% pos3 = 2*n_dens_std_C;
% neg3 = 2*n_dens_std_C;
%
% pos1 = pos*(pos1_0/pos1_0);
% neg1 = neg*(neg1_0/neg1_0);
% pos2 = pos*(pos2/pos1_0);
% neg2 = neg*(neg2/neg1_0);
% pos3 = pos*(pos3/pos1_0);
% neg3 = neg*(neg3/neg1_0);
%
% for i = 1:length(neg1)
%     if neg1(i) > n_dens_mean(i)
%         neg3(i) = (neg3(i)/neg1(i))*n_dens_mean(i);
%         neg2(i) = (neg2(i)/neg1(i))*n_dens_mean(i);
%         neg1(i) = n_dens_mean(i);
%     end
% end
%
% figure
% scatter([1,2,3,4,5,6,7,8,9,10],n_dens,'r*')
% hold on
% errorbar([1,2,3,4,5,6,7,8,9,10],n_dens_mean,neg1,pos1,'o','LineWidth',1)
% errorbar([1,2,3,4,5,6,7,8,9,10],n_dens_mean,neg2,pos2,'o','LineWidth',1)
% errorbar([1,2,3,4,5,6,7,8,9,10],n_dens_mean,neg3,pos3,'o','LineWidth',1)
% scatter([1,2,3,4,5,6,7,8,9,10],n_dens_mean,'bo')
% legend("Input","Plasma Depth","Einstein A","Calibration","State Density",'FontSize',18,'Location','northwest')
% set(gca,'xtick',[1:10],'xticklabel',states,'FontSize',20)
% xlabel('Paschen State Notation','FontSize',22)
% ylabel("Population Density (#/m^3)",'FontSize',22)
% title('Uncertainty Sources','FontSize',22)
% xlim([0.9 10.1])
% ylim([0 1.1*max(n_dens_mean+pos)])

% figure
% scatter(E,MC_mean./g')
% xlabel("Energy (eV)")
% ylabel("n/g (#/m^3)")
% title("Degeneracy Weighted Population Densities")

% T_BB = 1200+273;
% for i=1:length(lambda)
% BB_theoretical_intensity(i)= ((2*h*c^2)/((lambda(i)*10^-9)^5))*(1/(exp(h*c/((lambda(i)*10^-9)*k_B*T_BB))-1));
% end
% T_BB = 1200+273-3;
% for i=1:length(lambda)
% BB_theoretical_intensity_low(i)= ((2*h*c^2)/((lambda(i)*10^-9)^5))*(1/(exp(h*c/((lambda(i)*10^-9)*k_B*T_BB))-1));
% end
% T_BB = 1200+273+3;
% for i=1:length(lambda)
% BB_theoretical_intensity_high(i)= ((2*h*c^2)/((lambda(i)*10^-9)^5))*(1/(exp(h*c/((lambda(i)*10^-9)*k_B*T_BB))-1));
% end
% figure
% hold on
% plot(lambda, BB_theoretical_intensity_low)
% plot(lambda,BB_theoretical_intensity)
% plot(lambda,BB_theoretical_intensity_high)
% legend('T=1470K','T=1473K','T=1476K')
% figure
% plot(lambda, 100*abs(BB_theoretical_intensity_high-BB_theoretical_intensity_low)./BB_theoretical_intensity_high)
% title("% Diferrence between 1476K and 1470K Radiance")

% writematrix(n_dens_MC, '2022-09-15 - Spectra\MC Results\DoF\n_dens_all.txt')
% writematrix(n_dens_MC_A, '2022-09-15 - Spectra\MC Results\DoF\n_dens_A.txt')
% % writematrix(n_dens_MC_L, '2022-09-15 - Spectra\MC Results\12mL\n_dens_L.txt')
% writematrix(n_dens_MC_C, '2022-09-15 - Spectra\MC Results\DoF\n_dens_C.txt')

toc

%% Monte Carlo Function

function n_dens = MC_n_dens(lambda,Ar4p,NIST_lambda,A,A_accuracy,wdth,h,c,k_B,L,NIST_data,profile_type,t,A_px,eta_lambda,eta,dL,include,true_shape,assumed_shape)

if ismember("L",include)
    L_int = normrnd(L,dL);
else
    L_int = L;
end
L_int = L;
L_profile = linspace(-L_int/2,L_int/2,1000)';
if profile_type == "Uniform"
    profile = ones(1000,1);
elseif profile_type == "Gaussian"
    profile = exp(-L_profile.^2);
elseif profile_type == "Bimodal"
    loc = L*1/1.5;
    sigma = L/1.07;
    profile = exp(-((L_profile-loc)/sigma).^2)+exp(-((L_profile+loc)/sigma).^2);
    profile = profile/max(profile);
end

S_model = zeros(size(lambda));
p0 = [];
for i = 1:width(Ar4p)
    for j = 1:width(Ar4p.(i))
        NIST_ix = find(NIST_lambda==Ar4p.(i)(1,j));
        p0 = [p0, Ar4p.(i)(1,j)];
        A_n = A(NIST_ix);
        if ismember("A",include)
            A_n = normrnd(A(NIST_ix),A_accuracy(NIST_ix)*A(NIST_ix));
        end
        lambda_n = NIST_lambda(NIST_ix);
        n_j = Ar4p.(i)(2,j);
        if true_shape == "Gauss"
            G0 = [1,lambda_n,wdth];
            S_n_profile = (2.*sqrt(log(2))./(sqrt(pi).*wdth.*1E-9)).*(profile.*n_j).*(A_n./(4.*pi)).*((h.*c)./(lambda_n.*1E-9)).*Gauss(G0,lambda);
        elseif true_shape == "Lorentz"
            L0 = [1,lambda_n,wdth];
            S_n_profile = (2.*sqrt(log(2))./(sqrt(pi).*wdth.*1E-9)).*(profile.*n_j).*(A_n./(4.*pi)).*((h.*c)./(lambda_n.*1E-9)).*Lorentz(L0,lambda);
        elseif true_shape == "Voigt"
            V0 = [1,lambda_n,wdth/2,wdth/2];
            S_n_profile = (2.*sqrt(log(2))./(sqrt(pi).*wdth.*1E-9)).*(profile.*n_j).*(A_n./(4.*pi)).*((h.*c)./(lambda_n.*1E-9)).*Voigt(V0,lambda);
        end
        S_n = trapz(L_profile,S_n_profile,1);
        S_model = S_model + S_n;
    end
end

% Needs updating!!!!!!
d_r_meas = 2E-3;
x1 = 18.8E-2;
x2 = (6+20+14)*1E-2;
theta = 0*pi/180;
dTheta = 1*pi/180;
dr = sqrt(((x1/(1-theta)^2)*dTheta)^2+((1/(1-theta))*d_r_meas)^2+(-(x1*theta)/(theta-1)^2*dTheta)^2+((-theta/(1-theta))*d_r_meas)^2+(1*(3*d_r_meas))^2);
r_err = abs(normrnd(0,dr));
% dr_BB = 1E-3 + r_err;
% dr_GD = 1E-3 + r_err;
% r_BB = 62E-2;
% r_GD = 62E-2;
r = x1 + x2;
sigma_dOmega = 2*sqrt(2)*(r_err/r);
if ismember("C",include)
    err_dOmega = normrnd(1,sigma_dOmega);
else
    err_dOmega = 1;
end
S_model = S_model*err_dOmega;

T_BB = 1200+273;
a = -3;
b = 3;
if ismember("C",include)
    U_BB = (b-a).*rand(1) + a;
else
    U_BB = 0;
end
T_BB_MC = T_BB + U_BB;

S_model_BB = zeros(size(lambda));
for i=1:length(lambda)
    BB_theoretical_intensity= ((2*h*c^2)/((lambda(i)*10^-9)^5))*(1/(exp(h*c/((lambda(i)*10^-9)*k_B*T_BB))-1));
    BB_theoretical_intensity_MC = ((2*h*c^2)/((lambda(i)*10^-9)^5))*(1/(exp(h*c/((lambda(i)*10^-9)*k_B*T_BB_MC))-1));
    BB_err_MC = BB_theoretical_intensity_MC/BB_theoretical_intensity;
    S_model_BB(i) = S_model(i)*BB_err_MC;
end

lens_A = 25.4E-3;
lens_r = 58.8E-2;
Omega = ((lens_A/lens_r^2)/4*pi)/(1024*1024);

S_DC = zeros(size(lambda));
for i = 1:length(lambda)
    eta1 = interp1(eta_lambda,eta,lambda(i));
    S_DC(i) = ((2/eta1)*h*c)/((lambda(i)*1E-9)^2*A_px*Omega);
end

S_model_BB_DC = S_model_BB + S_DC;
e_camera = zeros(size(lambda));

for i = 1:length(lambda)
    eta1 = interp1(eta_lambda,eta,lambda(i));
    Phi = (S_model_BB_DC(i)*(lambda(i)*1e-9)^2*A_px*Omega)/(h*c);
    SNR = Phi*eta1*t/sqrt(Phi*0.12*t + 2*t + 40^2);
    e_camera(i) = normrnd(0,S_model_BB_DC(i)/SNR); % !!!!!!
end

S = S_model_BB_DC + e_camera;

%% Spectral Fitting

[pks,locs,wdths] = findpeaks(S,lambda,'MinPeakWidth',0.3);

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
if assumed_shape == "Gauss"
    initguess = [pks',locs',wdths'];
elseif assumed_shape == "Lorentz"
    initguess = [pks',locs',wdths'];
elseif assumed_shape == "Voigt"
    initguess = [pks',locs',wdths'/2,wdths'/2];
end

% Fit Gaussian to spectrum from initial guesses
[spectrum, fit_params, res] = SpectralFitting(lambda, S, initguess, assumed_shape);
n = length(pks); % number of transitions
a = fit_params(1:n); % transition peak height
p0 = fit_params(n+1:2*n); % center wavelength
w = fit_params(2*n+1:3*n); % width of line

% Calculate line intensities
for i = 1:length(pks)
    I(i) = (sqrt(pi)*a(i)*(w(i)*10^-9))/(2*sqrt(log(2)));
end
I(I<0) = 0;

% I: intensity of each transition line (area under curve)
% spectrum: fitted Gaussian result
% NIST_lambda: wavelengths of chosen transitions
% NIST_ix: indices of chosen transitions in NIST table
% A: Einstein A coefficients for chosen transitions

%% Calculations

% Calculate number densities
n_dens_trans = (4*pi*I.*(p0*10^-9))./(A.*h.*c.*L);

% n_dens_mean_all = [mean(n_dens_trans([22])),mean(n_dens_trans([16])),mean(n_dens_trans([14,19])),mean(n_dens_trans([10,15,21,24])),mean(n_dens_trans([9,13,23])),mean(n_dens_trans([8])),mean(n_dens_trans([4,12,20])),mean(n_dens_trans([3,6,18])),mean(n_dens_trans([2,5,11,17])),mean(n_dens_trans([1,7]))];
n_dens_mean_all = [mean(n_dens_trans([22])),mean(n_dens_trans([16])),mean(n_dens_trans([14,19])),mean(n_dens_trans([10,15,21])),mean(n_dens_trans([9,13])),mean(n_dens_trans([8])),mean(n_dens_trans([12,20])),mean(n_dens_trans([3,6,18])),mean(n_dens_trans([2,5,11,17])),mean(n_dens_trans([7]))];

n_dens_lumped = sum(n_dens_mean_all);

n_dens = [n_dens_trans,n_dens_lumped];

%

end