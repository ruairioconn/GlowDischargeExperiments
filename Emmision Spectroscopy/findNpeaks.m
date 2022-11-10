function [pks,locs,wdths] = findNpeaks(wavelength,intensity,N)
%findNpeaks Returns N largest peaks in dataset
%   [pks,locs,wdths] = findNpeaks(wavelength,intensity,N)

[pks, locs, wdths] = findpeaks(intensity, wavelength);

[~,I] = maxk(pks,N);

pks = pks(sort(I));
locs = locs(sort(I));
wdths = wdths(sort(I));

end