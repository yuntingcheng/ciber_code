function [nmap,rnmap,shotmap]=noise_realization(flight,inst,field,diff,g1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generate photon noise realization given 
%raw slope map and gain factor g1,g2.
%Ref: get_shotmap.m
%
%Input:
%   -flight: 36265/36277/40030
%   -inst: 1 or 2 (TM1, TM2)
%   -field number:(see get_fields)
%   -diff: 0 (full integration), 1 (1st-2nd diff)
%   -g1:gain factor g1 [(e-/s)/(ADU/fr)]
%Output:
%   -nmap:rnmap+shotmap (ADU/fr)
%   -rnmap:read noise realization map (ADU/fr)
%   -shotmap:sim photon noise map (ADU/fr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    slopepath=sprintf('%s/%d/%s/',...
            '/Users/ytcheng/ciber/data',flight,'slopedata');
    mypaths = get_paths(flight);
    fields=get_fields(flight,inst);
    cp=get_cal_params('flight',flight);
    frate=cp(inst).framerate;

    darkstatdir=strcat('/Users/ytcheng/ciber/doc/',...
            '20160906_NoiseRealization/darkstat/',num2str(flight),'/');
    load(strcat(darkstatdir,'TM',num2str(inst),'_darkps'),'darkps');

    %%% full integration case
    if diff==0
    load(sprintf('%sTM%d_%s.mat',...
        slopepath,inst,fields(field).name));
    shotmap=photonnoise_realization(rawmap,g1,fields(field).nfr,frate);
    rnmap=readnoise_realization(darkps(field).Clf2d_ave,7);
    nmap=rnmap+shotmap;
    
    %%% 1st 2nd halves difference    
    elseif diff==1
    load(sprintf('%sTM%d_%s_1st.mat',...
        slopepath,inst,fields(field).name));
    load(sprintf('%sTM%d_%s_2nd.mat',...
        slopepath,inst,fields(field).name));
    shot1=photonnoise_realization(rawmap1,g1,fields(field).nfrhalf,frate);
    shot2=photonnoise_realization(rawmap2,g1,fields(field).nfrhalf,frate);
    shotmap=(shot1-shot2)./sqrt(2);
    rnmap=readnoise_realization(darkps(field).Cld2d_ave,7);
    nmap=rnmap+shotmap;
    end
    
return 
