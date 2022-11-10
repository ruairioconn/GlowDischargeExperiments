function [yc,B,L] = BaselineCorrect(x,y)

[Cp,Sl,Ic] = ischange(y,'linear');                          % Detect Changes, Calculates Slopes (& Intercepts)
[Cts,Edg,Bin] = histcounts(Sl,1024);                        % Histogram Of Slopes
[Max,Binmax] = max(Cts);                                    % Find Largest Bin
LinearRegion = (Bin==Binmax);                               % Logical Vector Of Values Corresponding To Largest Number Of Slopes
B = polyfit(x(LinearRegion), y(LinearRegion), 1);           % Linear Fit
L = polyval(B, x);                                          % Evaluate
yc = y - L;                                                 % Detrend

end