function spec_array_out = spec_array_subtract(spec_array1,spec_array2)
%SPEC_ARRAY_SUBTRACT Subtracts two arrays of spectroscopy data
%spec_array_out = spec_array_subtract(spec_array1,spec_array2)
spec_array_out = spec_array1;
spec_array_out(1:length(spec_array1),4) = spec_array1(1:length(spec_array1),4) - spec_array2(1:length(spec_array2),4);
end

