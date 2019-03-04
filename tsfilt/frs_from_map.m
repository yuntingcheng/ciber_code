function frs=frs_from_map(slopemap,nfr,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%map a simulated time stream with a input slopemap.
%Input:
% slopemap: signal map [ADU/fr]
% nfr: number of frs
% addRN: add white read noise 
% sig_cds: sigma_CDS for read noise [e-/s](see Bock et al. 2013 table 2)
%           sig_cds=sqrt(2)sig_read(see Garnett & Forrest 1993)
%           use G1 in get_cal_params to convert to ADU/fr
% inst: 1--I band, 2--Hband
% offin: time stream offset
%
%Output:
% frs: simulated ts frames, size=nfr:1024:1024; [ADU]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
p.KeepUnmatched = true;

p.addRequired('slopemap');
p.addRequired('nfr');
p.addOptional('addRN',1,@isnumeric);
p.addOptional('sig_cds',10,@isnumeric);
p.addOptional('inst',1,@isnumeric);
p.addOptional('offin',zeros(size(slopemap)),@isnumeric); 
p.parse(slopemap,nfr,varargin{:});

slopemap = p.Results.slopemap;
nfr= p.Results.nfr;
addRN= p.Results.addRN;
sig_cds= p.Results.sig_cds;
inst= p.Results.inst;
offin= p.Results.offin;


if addRN
    if inst==1;sig_cds=10;else sig_cds=9;end
    cp=get_cal_params;
    sig_read=sig_cds/cp(inst).adu2e/sqrt(2);
end


frs=zeros(nfr,size(slopemap,1),size(slopemap,2));
for infr=1:nfr
    
    if addRN
        rnmap=normrnd(0,sig_read,size(slopemap));
        frs(infr,:,:)=slopemap.*infr+rnmap+offin;
    else
        frs(infr,:,:)=slopemap.*infr+offin;
    end
    
end

return