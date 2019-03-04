function [fit] = fit_beta_model(r_arr,prof_arr,err_arr,r1,r2)

model =  @(p,r) p(1).*(1 + (r./p(2)).^2).^(-3*p(3)/2);

sp = find(r_arr > r1 & r_arr < r2);
xData = r_arr(sp);
yData = prof_arr(sp);
dyData = err_arr(sp);
dxData = zeros(size(xData));

p = [1e-1, 10, 2];
op.LowerBound = [0 0 0.68];
op.UpperBound = [1 inf inf];
[fit.params,fit.dParams,fit.gof,fit.stddev] = ...
    fitChiSquare(xData,yData,model,p,dxData,dyData,op);
fit.fit_model = [model(fit.params,r_arr)];
fit.radius = [r1,r2];

return