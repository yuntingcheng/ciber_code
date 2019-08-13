function map = make_ihlmap_unisphere(flight,inst,Ifrac,m_min,m_max,rvir)
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

loaddir=strcat(mypaths.ciberdir,'doc/20170617_Stacking/psf_analytic/TM',...
    num2str(inst),'/');
load(strcat(loaddir,'fitpsfdat'),'fitpsfdat');

catdir=strcat(mypaths.ciberdir, 'doc/20170617_Stacking/maps/catcoord/SIDES/');
catfile=strcat(catdir,'sides.txt');

M = csvread(catfile,1);

x_arr=squeeze(M(:,2)');
y_arr=squeeze(M(:,3)');
x_arr=x_arr+1;
y_arr=y_arr+1;
rv_arr=squeeze(M(:,4)');

sr = ((7./3600.0)*(pi/180.0)).^2;
% use linear interpolated magnitude
if inst==1
    m_arr = squeeze(M(:,6)');
    lambdaeff=1.05;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
else
    m_arr = squeeze(M(:,7)');
    lambdaeff=1.79;
    I_arr=3631*10.^(-m_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
end

%%% select cat data %%%
sp=find(m_arr<=m_max & m_arr>m_min);
subI_arr=I_arr(sp);
subx_arr=x_arr(sp);
suby_arr=y_arr(sp);
subr_arr=rv_arr(sp).*rvir;

%%% get x,y coord in small grid %%%
xsmall_arr = subx_arr.*10 - 4.5;
ysmall_arr = suby_arr.*10 - 4.5;

%%% make srcmap w/o PSF %%%
Nsrc = numel(subI_arr);
srcmap=zeros(7200);
disp(sprintf('make IHLmap TM%d, mrange=(%d,%d), %d srcs',...
    inst,m_min,m_max,Nsrc));
if Nsrc>0
    Nlarge =7200+300+300;
    radmap = make_radius_map(zeros(2*Nlarge+1),Nlarge+1,Nlarge+1).*0.7;

    if Nsrc>20
        idx_print=floor(Nsrc/20) * (1:20);
    else
        idx_print=[];
    end

    print_count=0;
    for i=1:Nsrc
        xi = round(xsmall_arr(i));
        yi = round(ysmall_arr(i));
        rv = subr_arr(i);
        dx = Nlarge + 1 - xi;
        dy = Nlarge + 1 - yi;
        Imap = radmap(1+dx:7200+dx,1+dy:7200+dy);
        sp0 = find(Imap >= rv);
        sp1 = find(Imap < rv);
        Imap(sp0)=0;
        Imap(sp1)=sqrt(rv.^2 - Imap(sp1).^2);
        Imap(sp1)=subI_arr(i)*Ifrac.*Imap(sp1)./sum(Imap(sp1));
        srcmap = srcmap + Imap;

        if ismember(i,idx_print)
            print_count = print_count + 1;
            fprintf('make IHL of %d %% sources\n',print_count*5);
        end
    end
end
map = rebin_map_coarse(srcmap,10);
map = map.*100;
end