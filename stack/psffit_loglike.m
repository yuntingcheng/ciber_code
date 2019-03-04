function loglike=psffit_loglike(psfdata,m)
tmp = num2cell(m);
%[x0,y0,A,B,sig,r0,alpha] = deal(tmp{:});
[x0,y0,A,B] = deal(tmp{:});
sig=4.2;
r0=3.8;
alpha=3.6;
centx=ceil(size(psfdata,1)/2);
centy=ceil(size(psfdata,2)/2);
rcent_arr=make_radius_map(psfdata,centx,centy).*0.7;
sp=find(rcent_arr<20);
r_arr=make_radius_map(psfdata,centx-x0,centy-y0).*0.7;
psfmodel=A.*exp(-r_arr.^2./2./sig^2)+B./(1+(r_arr./r0).^alpha);
loglike=-(psfdata(sp)-psfmodel(sp)).^2./(r_arr(sp)+1)./psfmodel(sp).^2;
loglike=sum(loglike(:));
end