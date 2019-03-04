function [fit,fit_model] = fit_pl_prof_chi2(r_arr,prof_arr,err_arr,rmin, rmax)

sp = find(r_arr > rmin & r_arr<rmax);
xData = r_arr(sp);
yData = prof_arr(sp);
dyData = err_arr(sp);
dxData = zeros(size(xData));

% power law model function
model =  @(p,x) p(1).*x.^p(2);

% initial guess
p(1) = 1;
p(2) = -1.2;

op.LowerBound = [0 -inf];
op.UpperBound = [inf -1.01];
[fit.params,fit.dParams,fit.gof,fit.stddev] = ...
    fitChiSquare(xData,yData,model,p,dxData,dyData,op);

fit_model = model(fit.params,r_arr);
%fitmodel = model(fit.params, xData);
%errorbar(r_arr,prof_arr,err_arr,'b');hold on
%plot(xData,fitmodel,'r');
%set(gca,'XScale','log')

end