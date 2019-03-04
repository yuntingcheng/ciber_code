%{
date='03-13-2013';
band=1;
set_arr=1:5;
nfr_arr=2:23;
quad_corr=0;
%}
%%
%{
date='03-13-2013';
band=2;
set_arr=1:2;
nfr_arr=2:7;
quad_corr=1;
%}
%%
%{
date='04-25-2013';
band=2;
set_arr=1:3;
nfr_arr=2:6;
quad_corr=0;
%}
%%
%{
date='04-25-2013';
band=1;
set_arr=1:4;
nfr_arr=2:6;
quad_corr=0;
%}
%%
%{
date='03-22-2013';
band=2;
set_arr=1:4;
nfr_arr=2:6;
quad_corr=1;
%}
%%
%{
date='03-14-2013';
band=2;
set_arr=1:4;
nfr_arr=2:7;
quad_corr=1;
%}
%%
%{
date='05-13-2013';
band=1;
set_arr=1:3;
nfr_arr=2:5;
quad_corr=0;
%}
%%
%{
date='05-13-2013';
band=2;
set_arr=1:4;
nfr_arr=2:7;
quad_corr=0;
%}
%%
load(strcat('/Users/ytcheng/ciber/doc/20160808_DarkProcess/40030/band',...
    num2str(band),'_mask_inst'));
savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/Focus/fitTM'...
    ,num2str(band),'/'); 
for set=set_arr
    for nfr=nfr_arr
        framedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/',...
            'Focus/slopedata/',date,'/TM',num2str(band),...
                        '/set',num2str(set),'/');

        fitdat=fit_focus_g1(framedir,mask_inst,nfr,'quad_corr',quad_corr);
        if numel(fitdat)>0
            dataname=sprintf('%snfr%d_%s_set%d',savedir,nfr,date,set);
            save(dataname,'fitdat');
            disp(sprintf('%snfr%d_%s_set%d done',savedir,nfr,date,set));
        end
        
    end
end
