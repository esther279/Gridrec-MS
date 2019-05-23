function val = fun_legendre(coeff,x)

% the first coefficient is for 0-order

% legendreP(n,x) returns the nth degree Legendre polynomial at x.

val = 0;
for n = 0:length(coeff)-1
    val = val + coeff(n+1)*legendreP(n,x);
end

end

