clear;clc;

state_mapping = {[22];[16];[14,19];[10,15,21,24];[9,13,23];[8];[4,12,20];[3,6,18];[2,5,11,17];[1,7]};
Ar_2s4_transitions = [696.5431,706.7218,714.7042,763.5106,772.3761,801.4786,811.5311,912.2967];
Ar_2s3_transitions = [667.7282,727.2936,738.3980,751.4652,800.6157,810.3693,842.4648,965.7784];
Ar_2s2_transitions = [772.4207,794.8176,866.7944];
Ar_2s1_transitions = [750.3869,826.4522,840.8210,852.1442,922.4499,935.4220,978.5403];
Ar_4p_4s_transitions = sort([Ar_2s1_transitions,Ar_2s2_transitions,Ar_2s3_transitions,Ar_2s4_transitions]);
NIST_lambda = Ar_4p_4s_transitions(Ar_4p_4s_transitions>=650 & Ar_4p_4s_transitions<=970);

% density_input = readmatrix('2022-09-15 - Spectra\1 Torr densities.csv')*1E13;
% All_res = readmatrix('2022-09-15 - Spectra\MC Results\1 Torr\50W_uncertainties.txt');

density_input = readmatrix("C:\Users\Ruairi O'Connor\Box\Glow Discharge Data\2022-09-15 - Spectra\1 Torr densities.csv")*1E13;
All_res = readmatrix("C:\Users\Ruairi O'Connor\Box\Glow Discharge Data\2022-09-15 - Spectra\MC Results\1 Torr\50W_uncertainties.txt");

n_dens = density_input(10,:);
All_res(All_res<0) = 0;

MC_mean = mean(All_res,1);
% scatter_x = [1,2,3,3,4,4,4,4,5,5,5,6,7,7,7,8,8,8,9,9,9,9,10,10];
% scatter_y = MC_mean([22,16,14,19,10,15,21,24,9,13,23,8,4,12,20,3,6,18,2,5,11,17,1,7])*(100E-3/L_mean);

states = {"2p_1_0";"2p_9";"2p_8";"2p_7";"2p_6";"2p_5";"2p_4";"2p_3";"2p_2";"2p_1"};

% n_dens_mean_all = [mean(All_res([22])),mean(All_res([16])),mean(All_res([14,19])),mean(All_res([10,15,21,24])),mean(All_res([9,13,23])),mean(All_res([8])),mean(All_res([4,12,20])),mean(All_res([3,6,18])),mean(All_res([2,5,11,17])),mean(All_res([1,7]))];
n_dens_mean_all = [mean(All_res([22])),mean(All_res([16])),mean(All_res([14,19])),mean(All_res([10,15,21])),mean(All_res([9,13])),mean(All_res([8])),mean(All_res([12,20])),mean(All_res([3,6,18])),mean(All_res([2,5,11,17])),mean(All_res([7]))];

bins = linspace(floor(min(All_res,[],'all')),ceil(max(All_res,[],'all')),1000);
bins_plot = bins(1:end-1);

N = length(bins);
% hist_all = [histcounts(All_res(:,22),bins)/N;
%     histcounts(All_res(:,16),bins)/N;
%     histcounts(All_res(:,14),bins)/N + histcounts(All_res(:,19),bins)/N;
%     histcounts(All_res(:,10),bins)/N + histcounts(All_res(:,15),bins)/N + histcounts(All_res(:,21),bins)/N + histcounts(All_res(:,24),bins)/N;
%     histcounts(All_res(:,9),bins)/N + histcounts(All_res(:,13),bins)/N + histcounts(All_res(:,23),bins)/N;
%     histcounts(All_res(:,8),bins)/N;
%     histcounts(All_res(:,4),bins)/N + histcounts(All_res(:,12),bins)/N + histcounts(All_res(:,20),bins)/N;
%     histcounts(All_res(:,3),bins)/N + histcounts(All_res(:,6),bins)/N + histcounts(All_res(:,18),bins)/N;
%     histcounts(All_res(:,2),bins)/N + histcounts(All_res(:,5),bins)/N + histcounts(All_res(:,11),bins)/N + histcounts(All_res(:,17),bins)/N;
%     histcounts(All_res(:,1),bins)/N + histcounts(All_res(:,7),bins)/N;
%     histcounts(All_res(:,25),bins)/N];
hist_all = [histcounts(All_res(:,22),bins)/N;
    histcounts(All_res(:,16),bins)/N;
    histcounts(All_res(:,14),bins)/N + histcounts(All_res(:,19),bins)/N;
    histcounts(All_res(:,10),bins)/N + histcounts(All_res(:,15),bins)/N + histcounts(All_res(:,21),bins)/N;
    histcounts(All_res(:,9),bins)/N + histcounts(All_res(:,13),bins)/N;
    histcounts(All_res(:,8),bins)/N;
    histcounts(All_res(:,12),bins)/N + histcounts(All_res(:,20),bins)/N;
    histcounts(All_res(:,3),bins)/N + histcounts(All_res(:,6),bins)/N + histcounts(All_res(:,18),bins)/N;
    histcounts(All_res(:,2),bins)/N + histcounts(All_res(:,5),bins)/N + histcounts(All_res(:,11),bins)/N + histcounts(All_res(:,17),bins)/N;
    histcounts(All_res(:,7),bins)/N;
    histcounts(All_res(:,25),bins)/N];

for i=1:11
    hist_all(i,:) = hist_all(i,:)/trapz(bins(1:end-1),hist_all(i,:));
end

for i=1:11
    mu_state(i) = trapz(bins_plot,bins_plot.*hist_all(i,:));
    var_state(i) = trapz(bins_plot,(bins_plot-mu_state(i)).^2.*hist_all(i,:));
    sigma_state(i) = sqrt(var_state(i));
end

for i=1:11
    f = cumtrapz(bins_plot,hist_all(i,:));
    [~,loc] = min(abs(f-0.025));
    cent5(i) = bins_plot(loc);
    [~,loc] = min(abs(f-0.975));
    cent95(i) = bins_plot(loc);
end

hist_total = conv(hist_all(1,:),conv(hist_all(2,:),conv(hist_all(3,:),conv(hist_all(4,:),conv(hist_all(5,:),conv(hist_all(6,:),conv(hist_all(7,:),conv(hist_all(8,:),conv(hist_all(9,:),hist_all(10,:))))))))));
bins_conv = linspace(10*bins_plot(1),10*bins_plot(end),10*length(bins_plot)-9);
hist_total = hist_total/trapz(bins_conv,hist_total);

figure
t = tiledlayout(2,5);
title(t,'State PDFs','FontSize',22)
for i = 1:10
    nexttile
%     area(bins_plot(1:floor(end/15)),hist_all(i,1:floor(end/15)))
    area(bins_plot,hist_all(i,:))
    xline(mu_state(i),'r','LineWidth',2)
    xline(cent5(i),'LineWidth',2)
    xline(cent95(i),'LineWidth',2)
    xline(n_dens(i),'g','LineWidth',2)
    title(states(i),'FontSize',11)
    ylim([0 1.1*max(hist_all(i,1:floor(end/10)))])
    xlim([0 1E13])
    ax = gca; 
    ax.FontSize = 20; 
end
xlabel(t,'Population Density (#/m^3)','FontSize',22)
ylabel(t,'Probability Density','FontSize',22)

% figure
% scatter([1,2,3,4,5,6,7,8,9,10],n_dens,'r*')
% hold on
% dens_scatter = scatter(scatter_x,scatter_y);
% errorbar([1,2,3,4,5,6,7,8,9,10],mu_state,mu_state-cent5,cent95-mu_state,'bo','LineWidth',1)
% legend("Input","State Scatter","Monte Carlo Mean/StdDev",'FontSize',18,'Location','northwest')
% % legend("Monte Carlo Results",'FontSize',18,'Location','northwest')
% set(gca,'xtick',[1:10],'xticklabel',states,'FontSize',20)
% xlabel('Paschen State Notation','FontSize',22)
% ylabel("Population Density (#/m^3)",'FontSize',22)
% title('Monte Carlo Results','FontSize',22)
% xlim([0.9 10.1])

figure
area(bins_plot,hist_all(11,:))
xline(mu_state(11),'r','LineWidth',2)
xline(cent5(11),'LineWidth',2)
xline(cent95(11),'LineWidth',2)
xline(n_dens(11),'g','LineWidth',2)
title('Lumped State','FontSize',44)
xlabel('Population Density (#/m^3)','FontSize',44)
ylabel('Probability Density','FontSize',44)
xlim([0.25E13 7E13])
% ylim([0 4E-13])
ax = gca; 
ax.FontSize = 44; 

figure
% red blue dark green deep purple/magenta
newcolors = [0 0 1
    0 0.5 0
    1 0 0
    0.4940 0.1840 0.5560];
colororder(newcolors);
t = tiledlayout(2,5);
title(t,'State PDFs - Transitions','Fontsize',22)
for i=1:length(state_mapping)
    nexttile
    legend_str = {};
    plot_data = [];
    for j=1:length(state_mapping{i})
        plot(bins_plot,histcounts(All_res(:,state_mapping{i}(j)),bins),'LineWidth',2)
        hold on
        legend_str{j} = num2str(round(NIST_lambda(state_mapping{i}(j)),1)) + " nm";
        plot_data = [plot_data, histcounts(All_res(:,state_mapping{i}(j)),bins)];
    end
    title(states(i))
    legend(legend_str,'FontSize',20,'Location','northeast')
    ax = gca; 
    ax.FontSize = 20; 
%     xlim([0 bins_plot(floor(end/15))])
    xlim([0 1.5E13])
    ylim([0 1.25*max(plot_data)])
    legend boxoff
end
xlabel(t,'Population Density (#/m^3)','FontSize',22)
ylabel(t,'Probability Density','FontSize',22)
