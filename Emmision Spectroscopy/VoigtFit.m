function [fit,fitparams,res] = VoigtFit(x,y,initGuess)
%VoigtFit Finds best fit of spectrum to Voigt profile using least-squared
%optimization
%   INPUTS:
%       x: domain of fit
%       y: intensity data to fit
%       initGuess: initial guess for fit parameters
%   OUTPUTS:
%       fit: fitted spectrum
%       res: residual error

n = length(initGuess)/4;

for j = 1:n
    lb((j-1)*4 + 1) = initGuess((j-1)*4 + 1) - 0.5;
    lb((j-1)*4 + 2) = initGuess((j-1)*4 + 2)*0.9;
    lb((j-1)*4 + 3) = initGuess((j-1)*4 + 3)*0.9;
    lb((j-1)*4 + 4) = initGuess((j-1)*4 + 4)*0.95;
    
    ub((j-1)*4 + 1) = initGuess((j-1)*4 + 1) + 0.5;
    ub((j-1)*4 + 2) = initGuess((j-1)*4 + 2)*1.1;
    ub((j-1)*4 + 3) = initGuess((j-1)*4 + 3)*1.1;
    ub((j-1)*4 + 4) = initGuess((j-1)*4 + 4)*1.05;
end

fitparams = lsqcurvefit(@VoigtProfile, initGuess, x, y,lb,ub);
fit = VoigtProfile(fitparams,x);
res = fit - y;

end

