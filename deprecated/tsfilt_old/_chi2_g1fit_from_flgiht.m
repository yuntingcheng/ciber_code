function [chi2_arr,fCl,nCl_arr]=chi2_g1fit_from_flgiht...
    (rnmap_arr,flightmap,slopemap1,slopemap2,mask,nfr,frate,G1,varargin)

 %% Parse data
  p = inputParser;
  
  p.addRequired('rnmap_arr',@isnumeric);
  p.addRequired('flightmap',@isnumeric);
  p.addRequired('slopemap1',@isnumeric);
  p.addRequired('slopemap2',@isnumeric);
  p.addRequired('mask',@isnumeric);
  p.addRequired('nfr',@isnumeric);
  p.addRequired('frate',@isnumeric);
  p.addRequired('G1',@isnumeric);
  p.addOptional('pixscale',7,@isnumeric);
  p.addOptional('w',ones(1024),@isnumeric);
  p.addOptional('makeplot',0,@isnumeric);
  p.addOptional('saveplot',0,@isnumeric);
  p.addOptional('plotname','',@isstr);
  
  p.parse(rnmap_arr,flightmap,slopemap1,slopemap2,mask,nfr,frate,G1,varargin{:});

  rnmap_arr = p.Results.rnmap_arr;
  flightmap = p.Results.flightmap;
  slopemap1 = p.Results.slopemap1;
  slopemap2 = p.Results.slopemap2;
  mask = p.Results.mask;
  nfr = p.Results.nfr;
  frate = p.Results.frate;
  G1 = p.Results.G1;
  pixscale = p.Results.pixscale;
  w = p.Results.w;
  makeplot = p.Results.makeplot;
  saveplot = p.Results.saveplot;
  plotname = p.Results.plotname;
  clear p varargin;
%%
Nrn=size(rnmap_arr,1);
nCl_arr=zeros(Nrn,29);
rnCl_arr=zeros(Nrn,29);
phCl_arr=zeros(Nrn,29);
for i=1:Nrn
    rnmap=squeeze(rnmap_arr(i,:,:));
    phmap1=photonnoise_realization(slopemap1,G1,nfr,frate);
    phmap2=photonnoise_realization(slopemap2,G1,nfr,frate);
    phmap=(phmap1-phmap2)./2;
    nmap=(rnmap+phmap).*mask;
    nmap=nmap-mean(nmap(find(nmap)));nmap=nmap.*mask;
    [nCl] = get_angular_spec(nmap,nmap,pixscale,'w',w);
    nCl_arr(i,:)=nCl;
    
    rnmap=rnmap.*mask;
    rnmap=rnmap-mean(rnmap(find(rnmap)));rnmap=rnmap.*mask;
    [rnCl] = get_angular_spec(rnmap,rnmap,pixscale,'w',w);
    rnCl_arr(i,:)=rnCl;

    phmap=phmap.*mask;
    phmap=phmap-mean(phmap(find(phmap)));phmap=phmap.*mask;
    [phCl] = get_angular_spec(phmap,phmap,pixscale,'w',w);
    phCl_arr(i,:)=phCl;

end
avgnCl=squeeze(mean(nCl_arr));
varnCl=squeeze(var(nCl_arr));
avgrnCl=squeeze(mean(rnCl_arr));
avgphCl=squeeze(mean(phCl_arr));

fmap=flightmap.*mask;
fmap=fmap-mean(fmap(find(fmap)));fmap=fmap.*mask;
[fCl,l] = get_angular_spec(fmap,fmap,pixscale,'w',w);

chi2_arr=(fCl-avgnCl).^2./varnCl;
%%
if makeplot
    fig = figure;
    errorbar(l,l.*(l+1).*avgnCl./2./pi,...
        l.*(l+1).*sqrt(varnCl)./2./pi,'k.');hold on
    ax = get(fig,'CurrentAxes');
    set(ax,'XScale','log','YScale','log')
    loglog(l,l.*(l+1).*avgrnCl./2./pi,'b');
    loglog(l,l.*(l+1).*avgphCl./2./pi,'g');
    loglog(l,l.*(l+1).*fCl./2./pi,'r+');
    xlim([1e2,2e5]);
    legend({'RN+ph','RN','ph','flight'},'location','southeast');
    legend boxoff;
    title(sprintf('G1=%.2f',G1));
    drawnow
    if saveplot
    print(plotname,'-dpng');close
    end

end

end