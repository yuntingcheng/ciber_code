function  get_darkstat(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - make the 2DCl darkstat(full and diff) with only sigclip mask. 
% This 2DCl template is used for generate read noise in sims.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixscale=7;
mypaths=get_paths(flight);

loaddir=sprintf('%sTM%d/',mypaths.filtmap,inst);
load(sprintf('%s/darklongdat',loaddir),'darklongdat');
DCtemplate=darklongdat.DCtemplate_nomask; clear darklongdat

DCtempmask=zeros(1024);DCtempmask(find(DCtemplate))=1;
savedir=strcat(mypaths.alldat,'TM',num2str(inst),'/');
%% get and save dark stat
for ifield=4:8
disp(sprintf('ifield=%d',ifield));
dt=get_dark_times(flight,inst,ifield);
loaddir=sprintf('%sTM%d/field%d/',mypaths.filtmap,inst,ifield);

Cl2d_stack=zeros(numel(dt.time),1024,1024);
for i=1:numel(dt.time)
    load(strcat(loaddir,'labmap',num2str(i)),'labmap');
    
    filtmapf=labmap.filtmapf-DCtemplate;
    [~,maskf]=get_skymap(filtmapf,DCtempmask,4,5);
    filtmapf=filtmapf-mean(filtmapf(find(maskf)));filtmapf=filtmapf.*maskf;
    filtmapf=dc_offset_remove(filtmapf,maskf).*maskf;
    [~,~,~,~,~,~,Cl2d]=get_angular_spec(filtmapf,filtmapf,pixscale);
    Cl2d_stack(i,:,:)=Cl2d;
end
fCl2d_avg=squeeze(mean(Cl2d_stack));
fCl2d_std=squeeze(std(Cl2d_stack));

darkstat(ifield).Cl2d_avg=fCl2d_avg;
darkstat(ifield).Cl2d_std=fCl2d_std;

end
save(strcat(savedir,'darkstat'),'darkstat');

return