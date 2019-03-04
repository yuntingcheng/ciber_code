function [Cltot,l,lbins,Clsig,Clshot]=sigCl_extrap_mz14
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulated signal fitted from Fig S13 'output power spectrum' by eye
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,l,~,~,lbins]=get_angular_spec(randn(1024),randn(1024),7);
logx=[log10(2e3) log10(1e2)];x=10.^logx;
logy=[log10(2) log10(12)];y=10.^logy;
logCl=interp1(logx,logy,log10(l),'linear','extrap');
Clshot=ones(size(logCl)).*1e3.*2.*pi./1e5./(1e5+1);
Clsig=(((10.^logCl)).*2.*pi./l./(l+1));
Cltot = Clsig + Clshot;
return