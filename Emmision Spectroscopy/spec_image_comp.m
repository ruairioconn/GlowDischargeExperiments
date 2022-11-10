function comp_image = spec_image_comp(filedir)
%SPEC_IMAG_PROC converts series of spectroscopy images into composite array
%   Looks in project folder (determined by date) for images of certain
%   filetype (lamp, background, blackbody) and composites n images into
%   one. Then converts to an array

filePattern = fullfile(filedir, '*.csv'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
%     fprintf(1, 'Now reading %s\n', fullFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    if k == 1
        comp_image = readmatrix(fullFileName);
    else
        temp_image = readmatrix(fullFileName);
        comp_image(1:length(comp_image),4) = comp_image(1:length(comp_image),4) + temp_image(1:length(comp_image),4);
    end
end
comp_image(1:length(comp_image),4) = comp_image(1:length(comp_image),4)/length(theFiles);
end



