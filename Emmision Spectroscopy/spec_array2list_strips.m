function [wavelength, intensity] = spec_array2list_strips(spec_array)
%SPEC_ARRAY2LIST Converts array of spectroscopy data into wavelength and
%intensity lists
for i=1:length(spec_array)
    I_array(spec_array(i,2)+1,spec_array(i,1)+1)=spec_array(i,4);
    W_array(spec_array(i,2)+1,spec_array(i,1)+1)=spec_array(i,3);
end

h = height(I_array);
interval = floor(linspace(1,h,6));
for i = 1:length(interval)-1
    intensity(i,:) = mean(I_array(interval(i):interval(i+1)-1,:),1);
    wavelength(i,:) = mean(W_array(interval(i):interval(i+1)-1,:),1);
end

end