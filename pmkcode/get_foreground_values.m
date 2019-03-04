function foredat = get_foreground_values(flight,field)


%% Neutral Hydrogen (LABwavg)
% hydrogen col densities in cm^-2 from
% http://heasarc.nasa.gov/cgi-bin/Tools/w3nh/w3nh.pl
% weighted average for a 4 degree cone from LAB   

%% 100um Dust (iras)
% from an email from Matsuura on Aug 1st 2013,
% for the LRS FOV in MJy/sr

%% 25um DIRBE
% from Louis Levenson on Aug 15th 2013 in MJy/sr
%% to convert Wrights values from MJy/sr to nW/m2/sr
lam_m = 1.25*1e-6;
nu = 2.99e8/lam_m;
conv =1e9*1e6*1e-26*nu;

%% to convert Louis' 25um DIRBE data to nW/m2/sr
lam_m = 25*1e-6;
nu = 2.99e8/lam_m;
conv2 =1e9*1e6*1e-26*nu;

%% 36265

foregnd(1).field{1} = 'SWIRE';
foregnd(1).kel1p25(1) = 244.29;
foregnd(1).wright1p25(1) = .111*conv;
foregnd(1).LABwavg(1) = 1.11E+020;
foregnd(1).J11_ISL(1)= 30.3;
foregnd(1).iras(1)= 0.62;
foregnd(1).dirbe25(1)= 24.1*conv2;

foregnd(1).field{end +1} = 'NEP';
foregnd(1).kel1p25(end +1) = 235.51;
foregnd(1).wright1p25(end +1) = .107*conv;
foregnd(1).LABwavg (end +1)= 3.90E+020;
foregnd(1).J11_ISL(end +1)= 77.6;
foregnd(1).iras(end +1)= 2.50;
foregnd(1).dirbe25(end + 1)= 24.3*conv2;

foregnd(1).field{end +1} = 'elat10';
foregnd(1).kel1p25(end +1) = 487.32;
foregnd(1).wright1p25(end +1) = .222*conv;
foregnd(1).LABwavg (end +1)= 5.39E+020;
foregnd(1).J11_ISL(end +1)= -1;
foregnd(1).iras(end +1)= 4.17;
foregnd(1).dirbe25(end + 1)= 39.7*conv2;

foregnd(1).field{end +1} = 'elat30';
foregnd(1).kel1p25(end +1) = 363.82;
foregnd(1).wright1p25(end +1) = .168*conv;
foregnd(1).LABwavg (end +1)= 2.32E+020;
foregnd(1).J11_ISL(end +1)= -1;
foregnd(1).iras(end +1)= 1.65;
foregnd(1).dirbe25(end + 1)= 39.7*conv2;

foregnd(1).field{end +1} = 'BootesA';
foregnd(1).kel1p25(end +1) = 319.10;
foregnd(1).wright1p25(end +1) = .146*conv;
foregnd(1).LABwavg (end +1)= 1.09E+020;
foregnd(1).J11_ISL(end +1)= 14.7;
foregnd(1).iras(end +1)= 0.57;
foregnd(1).dirbe25(end + 1)= 25.1*conv2;

foregnd(1).field{end +1} = 'BootesB';
foregnd(1).kel1p25(end +1) = 327.93;
foregnd(1).wright1p25(end +1) = .150*conv;
foregnd(1).LABwavg (end +1)= 1.05E+020;
foregnd(1).J11_ISL(end +1)= 12.8;
foregnd(1).iras(end + 1)= 0.57;
foregnd(1).dirbe25(end + 1)= 27.9*conv2;

%% 36277

foregnd(2).field{1} = 'Lockman';
foregnd(2).kel1p25(1) = 330.23;
foregnd(2).wright1p25(1) = .155*conv;
foregnd(2).J11_ISL(1)= 15.5;
foregnd(2).LABwavg (end +1)= 6.12E+019;
foregnd(2).iras(1)= 0.39;
foregnd(2).dirbe25(1)= 34.8*conv2;

foregnd(2).field{end +1} = 'SWIRE';
foregnd(2).kel1p25(end +1) = 269.12;
foregnd(2).wright1p25(end +1) = .125*conv;
foregnd(2).J11_ISL(end +1)= 30.3;
foregnd(2).LABwavg(end +1) = 1.11E+020;
foregnd(2).iras(end + 1)= 0.61;
foregnd(2).dirbe25(end + 1)= 23.5*conv2;

foregnd(2).field{end +1} = 'NEP';
foregnd(2).kel1p25(end +1) = 280.55;
foregnd(2).wright1p25(end +1) = .129*conv;
foregnd(2).J11_ISL(end +1)= 77.6;
foregnd(2).LABwavg (end +1)= 3.90E+020;
foregnd(2).iras(end + 1)= 2.48;
foregnd(2).dirbe25(end + 1)= 23.5*conv2;

foregnd(2).field{end +1} = 'elat30';
foregnd(2).kel1p25(end +1) = 403.53;
foregnd(2).wright1p25(end +1) = .186*conv;
foregnd(2).J11_ISL(end +1)= 38.2;
foregnd(2).LABwavg (end +1)= 3.36E+020;
foregnd(2).iras(end + 1)= 2.63;
foregnd(2).dirbe25(end + 1)= 32.6*conv2;

foregnd(2).field{end +1} = 'BootesB';
foregnd(2).kel1p25(end +1) = 328.47;
foregnd(2).wright1p25(end +1) = .153*conv;
foregnd(2).J11_ISL(end +1)= 12.8;
foregnd(2).LABwavg (end +1)= 1.05E+020;
foregnd(2).iras(end + 1)= 0.65;
foregnd(2).dirbe25(end + 1)= -1;

foregnd(2).field{end +1} = 'BootesA';
foregnd(2).kel1p25(end +1) = 321.20;
foregnd(2).wright1p25(end +1) = .150*conv;
foregnd(2).J11_ISL(end +1)= 14.7;
foregnd(2).LABwavg (end +1)= 1.09E+020;
foregnd(2).iras(end + 1)= 0.65;
foregnd(2).dirbe25(end + 1)= -1;

%% 40030

foregnd(3).field{1} = 'DGL';
foregnd(3).kel1p25(1) = 255.44;
foregnd(3).wright1p25(1) = .116*conv;
foregnd(3).J11_ISL(1)= 43.2;
foregnd(3).LABwavg(1)= 4.52E+020;
foregnd(3).iras(1)= 2.41;
foregnd(3).dirbe25(1)= 33.36*conv2;

foregnd(3).field{end +1} = 'NEP';
foregnd(3).kel1p25(end +1) = 248.99;
foregnd(3).wright1p25(end +1) = .113*conv;
foregnd(3).J11_ISL(end +1)= 77.6;
foregnd(3).LABwavg(end +1)= 3.90E+020;
foregnd(3).iras(end + 1)= 2.51;
foregnd(3).dirbe25(end+1)= 26.0*conv2;

foregnd(3).field{end +1} = 'Lockman';
foregnd(3).kel1p25(end +1) = 436.13;
foregnd(3).wright1p25(end +1) = .198*conv;
foregnd(3).J11_ISL(end +1)= 15.5;
foregnd(3).LABwavg (end +1)= 6.12E+019;
foregnd(3).iras(end + 1)= 0.39;
foregnd(3).dirbe25(end+1)= -1;

foregnd(3).field{end +1} = 'elat10';
foregnd(3).kel1p25(end +1) = 561.54;
foregnd(3).wright1p25(end +1) = .253*conv;
foregnd(3).J11_ISL(end +1)= 11.7;
foregnd(3).LABwavg (end +1)= 1.74E+020;
foregnd(3).iras(end + 1)= 1.18;
foregnd(3).dirbe25(end+1)= 35.6*conv2;

foregnd(3).field{end +1} = 'elat30';
foregnd(3).kel1p25(end +1) = 404.19;
foregnd(3).wright1p25(end +1) = .187*conv;
foregnd(3).J11_ISL(end +1)= 9.34;
foregnd(3).LABwavg (end +1)= 0.90E+020;
foregnd(3).iras(end + 1)= 0.61;
foregnd(3).dirbe25(end+1)= 35.6*conv2;

foregnd(3).field{end +1} = 'BootesB';
foregnd(3).kel1p25(end +1) = 300.36;
foregnd(3).wright1p25(end +1) = .141*conv;
foregnd(3).J11_ISL(end +1)= 12.8;
foregnd(3).LABwavg (end +1)= 1.05E+020;
foregnd(3).iras(end + 1)= 0.67;
foregnd(3).dirbe25(end+1)= 31.1*conv2;

foregnd(3).field{end +1} = 'BootesA';
foregnd(3).kel1p25(end +1) = 294.33;
foregnd(3).wright1p25(end +1) = .138*conv;
foregnd(3).J11_ISL(end +1)= 14.7;
foregnd(3).LABwavg (end +1)= 1.09E+020;
foregnd(3).iras(end + 1)= 0.67;
foregnd(3).dirbe25(end+1)= 30.4*conv2;

foregnd(3).field{end +1} = 'SWIRE';
foregnd(3).kel1p25(end +1) = 244.39;
foregnd(3).wright1p25(end +1) = .113*conv;
foregnd(3).J11_ISL(end +1)= 30.3;
foregnd(3).LABwavg(end +1) = 1.11E+020;
foregnd(3).iras(end + 1)= 0.63;
foregnd(3).dirbe25(end+1)= 24.6*conv2;

%% get some indices to find your field

switch flight
    case 36265
        ind = 1;
    case 36277
        ind = 2;
    case 40030
        ind = 3;
end
    
%% string matching to pluck out what you want
nfields = numel(foregnd(ind).kel1p25);
fieldmatch = (1:nfields)*0;

for i=1:nfields
    fieldmatch(i) = strcmp(field,foregnd(ind).field(i));
end

%% make the returnable structure
spot = find(fieldmatch);
foredat.kel1p25 = foregnd(ind).kel1p25(spot);
foredat.wright1p25 = foregnd(ind).wright1p25(spot);
foredat.J11_ISL = foregnd(ind).J11_ISL(spot);
foredat.LABwavg = foregnd(ind).LABwavg(spot);
[DGLspec foredat.nbsdgl] = makeDGL((800:5000),foredat.LABwavg);
foredat.iras = foregnd(ind).iras(spot);
foredat.dirbe25 = foregnd(ind).dirbe25(spot);


%%

return