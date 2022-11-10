function [wavelength, intensity] = spec_array2list(spec_array)
%SPEC_ARRAY2LIST Converts array of spectroscopy data into wavelength and
%intensity lists
for i=1:length(spec_array)
    I_array(spec_array(i,2)+1,spec_array(i,1)+1)=spec_array(i,4);
    W_array(spec_array(i,2)+1,spec_array(i,1)+1)=spec_array(i,3);
end

intensity = mean(I_array,1);
wavelength = mean(W_array,1);

end

