flight=40030;
inst=1;
ifield=8;
quad='A';
dt=get_dark_times(flight,inst,ifield);

catdir='/Users/ytcheng/ciber/doc/20171018_stackihl/trilegal/';
catfile=strcat(catdir,dt.name,'_',quad,'_srcmapcoord.txt');
M = csvread(catfile);
m_arr=M(:,1);
xrand_arr=M(:,2);
yrand_arr=M(:,3);
navg=numel(m_arr);
%%
pixscale=7;
[~,~,~,~,lbin]=get_angular_spec(randn(1024),randn(1024),pixscale);
Clin=sigCl_extrap_mz14;

Nsub=80;
Ngrid=5120/Nsub;
sigmap = map_from_power_spec(lbin,Clin,Ngrid,Ngrid,Nsub*7*2,1);

sigmap=sigmap./std(sigmap(:));
sigmap(find(sigmap<-1))=-1;
sigmap=sigmap+1;

sigmap=sigmap./sum(sigmap(:))*navg;
nintmap=round(sigmap);


diff=sum(nintmap(:))-navg;
comp=randi([1,Ngrid^2],1,-diff);
addmap=zeros(Ngrid);
addmap(comp)=1;

nintmap=nintmap+addmap;
imageclip(nintmap);

%%
x_arr=[];
y_arr=[];
for ix=1:Ngrid
    for iy=1:Ngrid
        nsrc=nintmap(ix,iy);
        if nsrc>0
            for in=1:nsrc
                x_arr=[x_arr,ix];
                y_arr=[y_arr,iy];
            end
        end
        
    end
end

x_arr=x_arr+normrnd(0,0.3,1,navg);
y_arr=y_arr+normrnd(0,0.3,1,navg);

x_arr=x_arr*Nsub-1;
y_arr=y_arr*Nsub-1;

%%
figure
plot(x_arr,y_arr,'.');
figure
plot(xrand_arr,yrand_arr,'.');
%%
M1=zeros(size(M));
M1(:,1)=m_arr;

sp=randperm(navg);
M1(:,2)=x_arr(sp)';

sp=randperm(navg);
M1(:,3)=y_arr(sp)';

%%
catfile1=strcat(catdir,dt.name,'_',quad,'_srcmapcoord_corr.txt');
dlmwrite(catfile1, M1, 'precision', '%.3f');
