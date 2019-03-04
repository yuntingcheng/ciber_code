function [fit,fit_model] = fit_pl2_prof_chi2(r_arr,prof_arr,err_arr,r1,r2,r3)

% 1st fit
model1 =  @(p,x) p(1).*x.^p(2);

sp1 = find(r_arr > r1 & r_arr < r2);
xData = r_arr(sp1);
yData = prof_arr(sp1);
dyData = err_arr(sp1);
dxData = zeros(size(xData));

p(1) = 1;
p(2) = -2.3;
op.LowerBound = [0 -inf];
op.UpperBound = [inf -1.01];
[fit1.params,fit1.dParams,fit1.gof,fit1.stddev] = ...
    fitChiSquare(xData,yData,model1,p,dxData,dyData,op);

% 2nd fit
y_r2 = model1(fit1.params,r2);
model2 =  @(p,x) y_r2.*(x/r2).^p;

sp2 = find(r_arr > r2 & r_arr < r3);
xData = r_arr(sp2);
yData = prof_arr(sp2);
dyData = err_arr(sp2);
dxData = zeros(size(xData));

p = -1.;
op.LowerBound = [-inf];
op.UpperBound = [inf];
[fit2.params,fit2.dParams,fit2.gof,fit2.stddev] = ...
    fitChiSquare(xData,yData,model2,p,dxData,dyData,op);

fit_model = [model1(fit1.params,r_arr(find(r_arr<r2))) ...
    model2(fit2.params,r_arr(find(r_arr>r2)))];
%fit_model(find(r_arr>r2)) = fit_model(find(r_arr>r2)) + ...
%    model(fit2.params,r_arr(find(r_arr>r2)));
fit.params = [fit1.params fit2.params];
fit.radius = [r1,r2,r3];
return