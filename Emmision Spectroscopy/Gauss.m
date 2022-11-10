function G = Gauss(G0, lambda)
%Gauss Computes a Gaussian profile for given set of input lines
%   G0: Fit parameters in form [pks1,pks2,pks3,...,pksN,lambda1,lambda2,lambda3,...lambdaN,wdths1,wdths2,wdths3,...wdthsN]
%       

n_lines = size(G0,1);
pks = G0(1:n_lines);
locs = G0(n_lines+1:2*n_lines);
wdths = G0(2*n_lines+1:3*n_lines);

G = zeros(size(lambda));

for i = 1:n_lines
   x = (lambda-locs(i))./(wdths(i)/2);
   G = G + pks(i).*exp(-log(2).*x.^2);
end

end

