%% State Trend Plotting

clear;clc;

filename = '2022-09-15 - Spectra\Number Densities 15092022.xlsx';
[~,sheet_names] = xlsfinfo(filename);
state_mapping = {[22];[16];[14,19];[10,15,21,24];[9,13,23];[8];[4,12,20];[3,6,18];[2,5,11,17];[1,7]};
state_mapping = {[22];[16];[14,19];[10,15,21];[9,13];[8];[12,20];[3,6,18];[2,5,11,17];[7]};
state_mapping = {[22];[16];[14,19];[15,21];[9,13];[8];[12,20];[3,6,18];[2,5,17];[7]};
Ar_2s4_transitions = [696.5431,706.7218,714.7042,763.5106,772.3761,801.4786,811.5311,912.2967];
Ar_2s3_transitions = [667.7282,727.2936,738.3980,751.4652,800.6157,810.3693,842.4648,965.7784];
Ar_2s2_transitions = [772.4207,794.8176,866.7944];
Ar_2s1_transitions = [750.3869,826.4522,840.8210,852.1442,922.4499,935.4220,978.5403];
Ar_4p_4s_transitions = sort([Ar_2s1_transitions,Ar_2s2_transitions,Ar_2s3_transitions,Ar_2s4_transitions]);
NIST_lambda = Ar_4p_4s_transitions(Ar_4p_4s_transitions>=650 & Ar_4p_4s_transitions<=970);
states = {"2p_1_0";"2p_9";"2p_8";"2p_7";"2p_6";"2p_5";"2p_4";"2p_3";"2p_2";"2p_1"};

for k=1:numel(sheet_names)
  data{k} = readtable(filename,"Sheet",sheet_names{k});
end

colors = [0 0 1
    0 0.447 0.741
    0.4660 0.6740 0.1880
    0 1 0
    0.9290 0.6940 0.1250
    0.8500 0.3250 0.0980
    0.6350 0.0780 0.1840
    1 0 0
    0.4940 0.1840 0.5560
    1 0 1];
linestyles = ["-","--","-.",":"];

for i=1:length(data)
    x = data{i}.(1);
    x = x(~isnan(x));
    legend_text = [];
    figure
    title(sheet_names{i},'FontSize',22)
    hold on
    for j=1:10 % [5 6 10]
        for k = 1:length(state_mapping{j})
            y = data{i}.(state_mapping{j}(k)+2);
            y = y(y>1E10);
            plot(x,y,"Color",colors(j,:),'LineStyle',linestyles(k))
            legend_text = [legend_text, states(j) + ": " + num2str(NIST_lambda(state_mapping{j}(k))) + "nm"];
            if k == 1
                text(x(1),y(1),states(j) + " \rightarrow",'Color', colors(j,:),'HorizontalAlignment','right','FontSize',22)
            end
        end
    end
    legend(legend_text,'FontSize',18)
    xlabel("Voltage (V)",'FontSize',22)
    ylabel("Population Density (#/m^3)",'FontSize',22)
    legend boxoff
    ax = gca; 
    ax.FontSize = 20; 
    xlim([x(1)-15 x(10)+35])
end