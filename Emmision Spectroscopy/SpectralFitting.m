function [fit,fit_params,res,A] = SpectralFitting(lambda,intensity,init_guess,profile)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

opts = optimset('Display','off');

if profile == "Lorentz"
    n_lines = size(init_guess,1);
    pks = init_guess(1:n_lines);
    locs = init_guess(n_lines+1:2*n_lines);
    wdths = init_guess(2*n_lines+1:3*n_lines);
    
    lb_peaks = pks.*0.8;
    lb_locs = locs - 0.5;
    lb_wdths = wdths.*0.1;
    
    ub_peaks = pks.*1.2;
    ub_locs = locs + 0.5;
    ub_wdths = wdths.*1.9;
    
    lb = [lb_peaks, lb_locs, lb_wdths];
    ub = [ub_peaks, ub_locs, ub_wdths];
    
    fit_params = lsqcurvefit(@Lorentz, init_guess, lambda, intensity,lb,ub,opts);
    fit = Lorentz(fit_params,lambda);
    res = fit - intensity;

    n_lines = size(fit_params,1);
    pks = fit_params(1:n_lines);
    locs = fit_params(n_lines+1:2*n_lines);
    wdths = fit_params(2*n_lines+1:3*n_lines);

    for i=1:n_lines
        lambda_area = linspace(locs(i)-3*wdths(i),locs(i)+3*wdths(i),1000);
        line = Lorentz([pks(i),locs(i),wdths(i)],lambda_area);
        A(i) = trapz(lambda_area*1E-9,line);
    end
end

if profile == "Gauss"
    n_lines = size(init_guess,1);
    pks = init_guess(1:n_lines);
    locs = init_guess(n_lines+1:2*n_lines);
    wdths = init_guess(2*n_lines+1:3*n_lines);
    
    lb_peaks = pks.*0.8;
    lb_locs = locs - 0.5;
    lb_wdths = wdths.*0.1;
    
    ub_peaks = pks.*1.2;
    ub_locs = locs + 0.5;
    ub_wdths = wdths.*1.9;
    
    lb = [lb_peaks, lb_locs, lb_wdths];
    ub = [ub_peaks, ub_locs, ub_wdths];
    
    fit_params = lsqcurvefit(@Gauss, init_guess, lambda, intensity,lb,ub,opts);
    fit = Gauss(fit_params,lambda);
    res = fit - intensity;

    n_lines = size(fit_params,1);
    pks = fit_params(1:n_lines);
    locs = fit_params(n_lines+1:2*n_lines);
    wdths = fit_params(2*n_lines+1:3*n_lines);

    for i=1:n_lines
        lambda_area = linspace(locs(i)-3*wdths(i),locs(i)+3*wdths(i),1000);
        line = Gauss([pks(i),locs(i),wdths(i)],lambda_area);
        A(i) = trapz(lambda_area*1E-9,line);
    end
end

if profile == "Voigt"
    n_lines = size(init_guess,1);
    pks = init_guess(1:n_lines);
    locs = init_guess(n_lines+1:2*n_lines);
    sigma = init_guess(2*n_lines+1:3*n_lines);
    gamma = init_guess(3*n_lines+1:4*n_lines);
    
    lb_peaks = pks.*0.8;
    lb_locs = locs - 0.5;
    lb_sigma = sigma.*0.1;
    lb_gamma = gamma.*0.1;
    
    ub_peaks = pks.*1.2;
    ub_locs = locs + 0.5;
    ub_sigma = sigma.*1.9;
    ub_gamma = gamma.*1.9;
    
    Gauss_guess = [pks,locs,sigma];
    lb_Gauss = [lb_peaks, lb_locs, lb_sigma];
    ub_Gauss = [ub_peaks, ub_locs, ub_sigma];
    
    Gauss_params = lsqcurvefit(@Gauss, Gauss_guess,lambda,intensity,lb_Gauss,ub_Gauss,opts);
    
    pks = Gauss_params(1:n_lines);
    locs = Gauss_params(n_lines+1:2*n_lines);
    sigma = Gauss_params(2*n_lines+1:3*n_lines);
    gamma = init_guess(3*n_lines+1:4*n_lines)*0;
    
    lb_peaks = pks.*0.8;
    lb_locs = locs - 0.5;
    lb_sigma = sigma.*0.1;
    lb_gamma = gamma.*0.1;
    
    ub_peaks = pks.*1.2;
    ub_locs = locs + 0.5;
    ub_sigma = sigma.*1.9;
    ub_gamma = gamma.*1.9;
    
    lb = [lb_peaks, lb_locs, lb_sigma, lb_gamma];
    ub = [ub_peaks, ub_locs, ub_sigma, ub_gamma];
    
    fit_params = lsqcurvefit(@Voigt, init_guess, lambda, intensity,lb,ub,opts);
    fit = Voigt(fit_params,lambda);
    res = fit - intensity;

    n_lines = size(fit_params,1);
    pks = fit_params(1:n_lines);
    locs = fit_params(n_lines+1:2*n_lines);
    sigma = fit_params(2*n_lines+1:3*n_lines);
    gamma = fit_params(3*n_lines+1:4*n_lines);

    for i=1:n_lines
        lambda_area = linspace(locs(i)-3*(sigma(i)+gamma(i)),locs(i)+3*(sigma(i)+gamma(i)),1000);
        line = Voigt([pks(i),locs(i),sigma(i),gamma(i)],lambda_area);
        A(i) = trapz(lambda_area*1E-9,line);
    end
end


end

