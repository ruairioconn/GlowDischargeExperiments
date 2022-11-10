function [frequencies, amplitudes, phases] = harmonicAnalysis(signal, f_fund, N_harm, Fs)
%harmonicAnalysis: calculates frequencies, amplitudes and phases of
%harmonics of input signal
%NOTE: MATLAB FFT is based in cosine function
%   [frequencies, amplitudes, phases] = harmonicAnalysis(signal, f_fund, N_harm, Fs)
%   Takes input signal as 1D array, outputs frequencies in Hz, amplitudes
%   in units of input signal and phase in radians
%signal - input signal as 1D array
%f_fund - theoretical fundamental frequency
%N_harm - number of harmonics of interest
%Fs - Sampling frequency
% Information taken from: https://stackoverflow.com/questions/33962554/finding-the-phase-of-each-harmonics-using-fft

N = length(signal); % Number of samples
df = Fs/N; % Frequency resolution
sep = (f_fund/df)*0.9; % Minimum separation when looking for harmonics

Xf = fft(signal); % take FFT
Nmax = N/2 + 1;
Xf = Xf(1:Nmax); % Make positive

ind = find(abs(Xf) ~= 0); % find indices where FFT ~= 0;
[~,idx] = findpeaks(abs(Xf(ind)),'MinPeakDistance',sep); %Find harmonics using peak finder with minimum separation distance
% Remove potential harmonics < fundamental frequency
if idx(1) < sep
    idx = idx(2:end);
end
% Set number of harmonics of interest
if length(idx) > N_harm
    idx = idx(1:N_harm);
end

Amp   = zeros(1,length(idx));
freq  = zeros(1,length(idx));
phase = zeros(1,length(idx));
for i=1:length(idx)
    ratio = abs(Xf(idx(i)+1))/abs(Xf(idx(i)));
    if (abs(Xf(idx(i)+1)) > abs(Xf(idx(i)-1)))
        ratio = -ratio;
    end
    
    frac = fractional_frequency(ratio, N);
    freq(i)  = (idx(i)-1+frac)*Fs/N;
    Amp(i)   = 2 * abs(Xf(idx(i))) * abs(sin(pi*frac/N)/sin(pi*frac));
    phase(i) = angle( Xf(idx(i)) .* (1-exp(2*pi*frac*1i/N)) ./ (1-exp(2*pi*frac*1i)) );
end

frequencies = freq;
amplitudes = Amp;
phases = phase;

% Solve for "f" for which ratio = sin(pi*frac/N)/sin(pi*(frac-1)/N)
    function f = fractional_frequency(ratio, N)
        
        niter = 20;
        K = (pi/N) * sin(pi/N);
        f = 0;
        for j=1:niter
            a  = sin(pi*f/N);
            b  = sin(pi*(f-1)/N);
            
            y  = ratio - a/b;
            yp = K / (b^2);
            f = max(-0.5, min(f - y/yp, 0.5));
        end
    end


end

