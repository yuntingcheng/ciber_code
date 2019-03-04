function map = make_srcmap_sides(flight,inst,ifield,m_min,m_max,finePSF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim src map SIDES
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - PSF: True-use full PSF, False-just put source in the center pix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,'sides.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

sr = ((7./3600.0)*(pi/180.0)).^2;
% use linear interpolated magnitude
if inst==1
    m_arr = squeeze(M(:,5)');
    lambdaeff=1.05;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
else
    m_arr = squeeze(M(:,6)');
    lambdaeff=1.79;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
end

%%% select cat data %%%
sp=find(m_arr<=m_max & m_arr>m_min);
subI_arr=I_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);

%%% get x,y coord in small grid %%%
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;

%%% get psf params %%%
beta=fitpsfdat(ifield).psfmodel.beta_best;
rc=fitpsfdat(ifield).psfmodel.rc_best;
norm=fitpsfdat(ifield).psfmodel.norm;

%%% make srcmap with PSF %%%
if finePSF
Nsrc = numel(subI_arr);
srcmap=zeros(7200);
disp(sprintf('make srcmap TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));
if Nsrc>0
    Nlarge =10240+300+300;
    radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;
    Imap_large = norm .* (1 + (radmap/rc).^2).^(-3.*beta./2);

    if Nsrc>20
        idx_print=floor(Nsrc/20) * (1:20);
    else
        idx_print=[];
    end

    print_count=0;
    for i=1:Nsrc
        xi = round(xsmall_arr(i));
        yi = round(ysmall_arr(i));
        dx = Nlarge + 1 - xi;
        dy = Nlarge + 1 - yi;
        Imap = Imap_large(1+dx:7200+dx,1+dy:7200+dy).*subI_arr(i);
        srcmap = srcmap + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            fprintf('stack %d %% sources\n',print_count*5);
        end
    end
end
map = rebin_map_coarse(srcmap,10);
map = map.*100;

%%% usePSF == Flase %%%
else
Nsrc = numel(subI_arr);
disp(sprintf('make srcmap TM%d %s, mrange=(%d,%d), %d srcs',...
    inst,dt.name,m_min,m_max,Nsrc));
map = zeros(720);
Nlarge =1024+30+30;
radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*7;
Imap_large = norm .* (1 + (radmap/rc).^2).^(-3.*beta./2);
Imap_large = Imap_large./sum(Imap_large(:));
if Nsrc>0
    if Nsrc>20
        idx_print=floor(Nsrc/20) * (1:20);
    else
        idx_print=[];
    end

    print_count=0;
    for i=1:Nsrc
        xi = round(subx_arr(i));
        yi = round(suby_arr(i));
        dx = Nlarge + 1 - xi;
        dy = Nlarge + 1 - yi;
        Imap = Imap_large(1+dx:720+dx,1+dy:720+dy).*subI_arr(i);
        map = map + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            fprintf('stack %d %% sources\n',print_count*5);
        end
    end
end

end

end