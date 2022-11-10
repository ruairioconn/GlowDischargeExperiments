function out = call_NIST(wavelengths, varargin)
%call_NIST Return parameters from NIST database for given wavelength range
%   [lambda, E, g, A, accuracy] = call_NIST(wavelengths, peak_locs)
%   Loads NIST url for wavelength range given, finds data for closest
%   matching peaks to given peak_locs


lambda_min = min(wavelengths);
lambda_max = max(wavelengths);

url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Ar+I&limits_type=0&low_w='+string(lambda_min)+'&upp_w='+string(lambda_max)+'&unit=1&submit=Retrieve+Data&de=0&format=3&line_out=0&remove_js=on&en_unit=1&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on';
% url = 'https://physics.nist.gov/cgi-bin/ASD/lines1.pl?spectra=Ar&limits_type=0&low_w='+string(lambda_min)+'&upp_w='+string(lambda_max)+'&unit=1&submit=Retrieve+Data&de=0&format=3&line_out=0&remove_js=on&en_unit=1&output=0&bibrefs=1&page_size=15&show_obs_wl=1&show_calc_wl=1&unc_out=1&order_out=0&max_low_enrg=&show_av=2&max_upp_enrg=&tsb_value=0&min_str=&A_out=0&intens_out=on&max_str=&allowed_out=1&forbid_out=1&min_accur=&min_intens=&conf_out=on&term_out=on&enrg_out=on&J_out=on&g_out=on';

% NIST = webread(url);
% NIST_str = convertCharsToStrings(NIST);
% fid = fopen('NIST_data_temp.txt','wt');
% fprintf(fid, NIST_str);
% fclose(fid);

NIST_accuracies = {'AAA','AA','A+','A','B+','B','C+','C','D+','D','E'};
NIST_accuracy_vals = [0.003,0.01,0.02,0.03,0.07,0.1,0.18,0.25,0.4,0.5,0.5];
accuracy_key = containers.Map(NIST_accuracies, NIST_accuracy_vals);

NIST_file = 'NIST_data_temp.txt';
opts = detectImportOptions(NIST_file);
NIST_data = readtable(NIST_file,opts);
header = {'lambda','ritz_lambda','intensity','A_ki','accuracy','E_i','E_k','conf_i','term_i','J_i','conf_k','term_k','J_k','g_i','g_k','Type','tp_ref','line_ref'};
% header = {'element','sp_num','lambda','lambda_uncertainty','ritz_lambda','ritz_lambda_uncertainty','intensity','A_ki','accuracy','E_i','E_k','conf_i','term_i','J_i','conf_k','term_k','J_k','g_i','g_k','Type','tp_ref','line_ref'};
NIST_data.Properties.VariableNames = header;

if isa(varargin{1},'numeric')
    peak_locs = varargin{1};
    for i = 1:length(peak_locs)
        [minValue, NIST_selection(i)] = min(abs(NIST_data.lambda-peak_locs(i)));
        if i > 1 && NIST_selection(i) == NIST_selection(i-1)
            oldIndex = NIST_selection(i);
            [minValue, NIST_selection(i)] = min(abs(NIST_data.lambda(oldIndex+1:end)-peak_locs(i)));
            NIST_selection(i) = NIST_selection(i) + oldIndex;
        end
    end

    N = length(NIST_selection);
    lambda = NIST_data.lambda(NIST_selection);
    E_upper = NIST_data.E_k(NIST_selection);
    E_lower = NIST_data.E_i(NIST_selection);
    g = NIST_data.g_k(NIST_selection);
    A = NIST_data.A_ki(NIST_selection);

    for i=1:length(NIST_selection)
        err_key = NIST_data.accuracy(NIST_selection(i));
        accuracy(i) = accuracy_key(err_key{1});
    end

    out.lambda = lambda;
    out.E_upper = E_upper;
    out.E_lower = E_lower;
    out.g = g;
    out.A = A;
    out.accuracy = accuracy;

else
    out = NIST_data;
end

end