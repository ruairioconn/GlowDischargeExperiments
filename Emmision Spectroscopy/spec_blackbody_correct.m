function [calibrated_intensity, BB_theoretical_intensity, scale_factor] = spec_blackbody_correct(data_wavelength,data_intensity,BB_intensity, temperature)
%SPEC_BLACKBODY_CORRECT
% [calibrated_intensity, BB_theoretical_intensity, scale_factor] = spec_blackbody_correct(data_wavelength,data_intensity,BB_intensity, temperature)
% [data_wavelength,data_intensity] = spec_array2list(data_array);
% [BB_wavelength,BB_intensity] = spec_array2list(BB_array);

h = 6.626E-34;
c = 299.792E6;
k_B = 1.38E-23;
T_BB = temperature+273;

for i=1:length(data_wavelength)
   BB_theoretical_intensity(i) = ((2*h*c^2)/((data_wavelength(i)*10^-9)^5))*(1/(exp(h*c/((data_wavelength(i)*10^-9)*k_B*T_BB))-1));
end

%BB_theoretical_intensity = BB_theoretical_intensity/max(BB_theoretical_intensity);

for i=1:length(BB_intensity)
   scale_factor(i) = BB_theoretical_intensity(i)/BB_intensity(i);
end

calibrated_intensity = scale_factor.*data_intensity;

end

