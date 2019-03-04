function map = make_srcmap_ps(flight,inst,ifield,type,m_min,m_max,interp)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the sim src map PanSTARR
%
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star, 0:all, 2:undefined
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mypaths=get_paths(flight);

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/TM',...
    num2str(inst),'/PanSTARRS/');
loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

dt=get_dark_times(flight,inst,ifield);

%%% read cat data %%%
catfile=strcat(catdir,dt.name,'.txt');

M = csvread(catfile,1);


x_arr=squeeze(M(:,4)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;

my_arr=squeeze(M(:,9)');
cls_arr=squeeze(M(:,10)');

sr = ((7./3600.0)*(pi/180.0)).^2;

% use y band mag
if interp == 0
    lambdaeff=0.9633;
    I_arr=3631*10.^(-my_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    m_arr = my_arr;
end

% use linear interpolated magnitude (but bin by y mag)
if interp == 1
    if inst==1
        mcb_arr = squeeze(M(:,11)');
        lambdaeff=1.05;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    else
        mcb_arr = squeeze(M(:,12)');
        lambdaeff=1.79;
        I_arr=3631*10.^(-mcb_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
    end
    m_arr = my_arr;
end

% use corrected magnitude
if interp ==2
    
    if inst==1
        mlin_arr = squeeze(M(:,11)');
    else
        mlin_arr = squeeze(M(:,12)');
    end
[m_arr, I_arr] = get_corrected_mag(inst, mlin_arr, my_arr, cls_arr);
end

%%% select cat data %%%
sp=find(m_arr<=m_max & m_arr>m_min & cls_arr==type ...
    & x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);

if type == 0
    sp=find(m_arr<=m_max & m_arr>m_min ...
         & x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);
end

if type == 2
    sp=find(m_arr<=m_max & m_arr>m_min & cls_arr~=1 & cls_arr~=-1 ...
         & x_arr>0.5-30 & x_arr<1024.5+30 & y_arr>0.5-30 & y_arr<1024.5+30);
end
subm_arr=m_arr(sp);
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

%%% make srcmap %%%
Nsrc = numel(subI_arr);
srcmap=zeros(10240);
disp(sprintf('make srcmap TM%d %s, mrange=(%d,%d), type = %d, %d srcs',...
    inst,dt.name,m_min,m_max,type,Nsrc));
if Nsrc>0
    Nlarge = 10240+300+300;
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
        Imap = Imap_large(1+dx:10240+dx,1+dy:10240+dy).*subI_arr(i);
        srcmap = srcmap + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            disp(sprintf('make srcmap %d %% sources',print_count*5));
        end
    end
end
map = rebin_map_coarse(srcmap,10);
map = map.*100;

end