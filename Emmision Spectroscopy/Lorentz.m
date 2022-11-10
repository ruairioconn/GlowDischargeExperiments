function L = Lorentz(L0, lambda)
%Lorentz Computes a Lorentz profile for given set of input lines
%   L0: Fit parameters in form [pks1,pks2,pks3,...,pksN,lambda1,lambda2,lambda3,...lambdaN,wdths1,wdths2,wdths3,...wdthsN]
%       

n_lines = size(L0,1);
pks = L0(1:n_lines);
locs = L0(n_lines+1:2*n_lines);
wdths = L0(2*n_lines+1:3*n_lines);

L = zeros(size(lambda));

for i = 1:n_lines
   x = (lambda-locs(i))/(wdths(i)/2);
   L = L + pks(i)./(1+x.^2);
end

end

