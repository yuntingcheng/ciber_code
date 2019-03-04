%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get PS flat bias factor for signal
% assuming signal is the same in all field, this can be 
% derived by hand
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
savedirfull='/Users/ytcheng/ciber/doc/20170209_TsFilter/fullPS/';
load(strcat(savedirfull,'flightCldat'),'flightCldat');

for ifield=[8,7,6,5,4,2,1]
    FFuse=[0,0,0,1,1,1,1,1];FFuse(ifield)=0;
    nusum=0;
    desum=0;
    for jfield=[8,7,6,5,4,2,1]
        if FFuse(jfield)==1
            nusum=nusum+(1/flightCldat(jfield).meanmap);
            desum=desum+sqrt(flightCldat(jfield).meanmap);
        end
    end
    errfac=flightCldat(ifield).meanmap^2*nusum/(desum^2);
    FFsigbias=1+errfac;
    flightCldat(ifield).FFsigbias=FFsigbias;
end

save(strcat(savedirfull,'flightCldat'),'flightCldat');