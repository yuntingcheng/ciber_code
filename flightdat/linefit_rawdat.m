for flight=[36265,36277,40030]
savedir=strcat('/Users/ytcheng/ciber/data/',num2str(flight),'/slopedata/');

for band=[1,2]

ft = get_field_times(flight,band);

switch flight
case 36265;name_arr=ft.name(14:4:end-4);name_arr(3)=[];
case 36277;name_arr=ft.name(1:3:end-2);
case 40030;name_arr=ft.name(2:3:end-3);
end


%%%%%%%%%%
for i=1:numel(name_arr)
name=name_arr{i};
[frames] = get_data_frames(band,name,'flight',flight,'verbose',0);

frames=frames(3:end,:,:);
nfr=size(frames);nfr=nfr(1);nfrhalf=floor(nfr/2);
pr=sprintf('%d,TM%d,%s,nfr=%d',flight,band,name,nfr);disp(pr);

[rawmap,~,rawmapvar]=linfit_map(frames);
[rawmap1,~,rawmap1var]=linfit_map(frames(1:nfrhalf,:,:));
[rawmap2,~,rawmap2var]=linfit_map(frames(1+nfrhalf:2*nfrhalf,:,:));

figure
subplot(2,2,1)
imageclip(rawmap);title(name)
subplot(2,2,2)
imageclip(rawmap1);title(name)
subplot(2,2,3)
imageclip(rawmap2);title(name)
subplot(2,2,4)
imageclip(rawmap1-rawmap2);title(name)

save(sprintf('%sTM%d_%s',savedir,band,name),'rawmap');
save(sprintf('%sTM%d_%s_1st',savedir,band,name),'rawmap1');
save(sprintf('%sTM%d_%s_2nd',savedir,band,name),'rawmap2');

save(sprintf('%sTM%d_%s_var',savedir,band,name),'rawmapvar');
save(sprintf('%sTM%d_%s_1st_var',savedir,band,name),'rawmap1var');
save(sprintf('%sTM%d_%s_2nd_var',savedir,band,name),'rawmap2var');
end

end
end
