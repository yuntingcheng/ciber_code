%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the NIST spectrum data and save in a struc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

band=1;

savedir=strcat('/Users/ytcheng/ciber/doc/20160912_CalFac/NIST/TM'...
                    ,num2str(band),'/');                

scandir=strcat('/Volumes/HD1TB/CIBER/data/Cal_NIST/TM'...
    ,num2str(band),'/NIST/');
scanfile=dir(strcat(scandir,'*.txt'));

for filenum=1:numel(scanfile)
filename=strcat(scandir,scanfile(filenum).name);

dat_arr = importdata(filename,' ');dat_arr=dat_arr(4:end,:);
Wavelength_arr=zeros(1,numel(dat_arr));
MeanSig_arr=zeros(1,numel(dat_arr));
GainCorrSig_arr=zeros(1,numel(dat_arr));
CalCoeff_arr=zeros(1,numel(dat_arr));
L_arr=zeros(1,numel(dat_arr));
for i=1:numel(dat_arr)
    dat=str2num(dat_arr{i});
    Wavelength_arr(i)=dat(1);MeanSig_arr(i)=dat(2);
    GainCorrSig_arr(i)=dat(3);CalCoeff_arr(i)=dat(4);L_arr(i)=dat(5);
end

NISTdat(filenum).name=scanfile(filenum).name;
NISTdat(filenum).Wavelength_arr=Wavelength_arr;
NISTdat(filenum).MeanSig_arr=MeanSig_arr;
NISTdat(filenum).GainCorrSig_arr=GainCorrSig_arr;
NISTdat(filenum).CalCoeff_arr=CalCoeff_arr;
NISTdat(filenum).L_arr=L_arr;
end
save(strcat(savedir,'NISTdat'),'NISTdat');
