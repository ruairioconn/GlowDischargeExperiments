function V = Voigt(V0, lambda)
%Voigt Computes a Voigt profile for given set of input lines
%   V0: Fit parameters in form [pks1,pks2,pks3,...,pksN,lambda1,lambda2,lambda3,...lambdaN,sigma1,sigma2,sigma3,...sigmaN,gamma1,gamma2,gamma3,...,gammaN]
%       

n_lines = size(V0,1);
pks = V0(1:n_lines);
locs = V0(n_lines+1:2*n_lines);
sigma = V0(2*n_lines+1:3*n_lines)/2;
gamma = V0(3*n_lines+1:4*n_lines)/2;

V = zeros(size(lambda));

for i = 1:n_lines
   x = lambda - locs(i);
   z = (x + 1i.*gamma(i))./(sigma(i).*sqrt(2));
   V = V + pks(i).*real(Faddeeva_w(z))./(sigma(i).*sqrt(2.*pi));
end

end

