function [L, r] = LorentzFit(wavelength, intensity, pks, locs, w)
%LORENTZFIT Function to calculate Lorentz fit to spectral data
n = length(wavelength);
m = length(pks);
c = pks;

L = zeros(1,length(wavelength));
for i=1:m
x = (wavelength-locs(i))./(w(i)/2);
L_i = 1./(1+x.^2);
L = L + c(i)*L_i;
end

for i=1:n
   r(i) = intensity(i) - L(i); 
end

min_function = 0;
for i=1:m
    x = (wavelength-locs(i))./(w(i)/2);
    min_function = min_function - 2.*r(i).*(1./(1+x.^2)+c(i).*(((2.*x)./(1+x.^2).^2).*((1)./(w(i)./2)))+c(i).*((4.*x.*(wavelength-locs(i)))./((1-x.^2).^2.*(w(i).^2))));
end
min_function;
end

