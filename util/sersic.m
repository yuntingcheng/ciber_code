function I_arr = sersic(x_arr,Ie,n,xe)
% https://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Sersic2D.html
% return Sersic profile: 
% I(r_arr) = Ie*exp(-bn((r_arr/re)^(1/n)-1));
% where bn is defined such that Gamma(2n) = 2 gamma(bn,2n);

% find bn, since bn~2n-1/3 (c.f. wiki) use this to set the test range
bn_arr = (2*n-(1/3)).*logspace(-1,1,1000);
if any(bn_arr<0)
    bn_arr = logspace(-30,1,1000);
end

% since gammainc(x,a) = gamma(x,a)/Gamma(a);
diff = gammainc(bn_arr, 2*n)-0.5;
sp = find(abs(diff) == min(abs(diff)));
bn = bn_arr(sp);

I_arr = (10^Ie) * exp(-bn*((x_arr/xe).^(1/n)-1));

return