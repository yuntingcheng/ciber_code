function plot_chi2g1fit(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the Chi2-G1 of each field, and chi2 fit with 
%all the fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if inst==1
    G1_FF=-2.68;
elseif inst==2
    G1_FF=-2.84;
end

fields_info=get_fields(flight,inst);
savedir=strcat('/Users/ytcheng/ciber/doc/20170325_alldat/TM',...
    num2str(inst),'/noisemodel/diffG1fit/');
load(strcat(savedir,'chi2dat'),'chi2dat');
G1_arr=chi2dat.G1_arr;
nfr_arr=chi2dat.nfr_arr;

figure
subchi2_allfield=zeros(size(G1_arr));
subdof_allfield=-1;

for ifield=4:8
dt=get_dark_times(flight,inst,ifield);
name=dt.name;

chi2totsum=zeros(1,numel(G1_arr));
subchi2totsum=zeros(1,numel(G1_arr));
doftot=-1;
dofsub=-1;

for nfr=nfr_arr
    field_use=chi2dat.nfr(nfr).field_use;
    chi2_arr=chi2dat.nfr(nfr).chi2tot_arr;
    subchi2_arr=chi2dat.nfr(nfr).subchi2tot_arr;
    if sum(find(field_use==ifield))
        ind=find(field_use==ifield);
        chi2=chi2_arr(ind,:);
        subchi2=subchi2_arr(ind,:);
        chi2totsum=chi2totsum+chi2;
        subchi2totsum=subchi2totsum+subchi2;
        doftot=doftot+21;
        dofsub=dofsub+5;
    end  
end
subchi2_allfield=subchi2_allfield+subchi2totsum;
subdof_allfield=subdof_allfield+dofsub;

plot(G1_arr,subchi2totsum/dofsub,'o-','color',get_color(ifield),...
    'DisplayName',name);hold on
minind=find(subchi2totsum==min(subchi2totsum));
plot(G1_arr(minind),subchi2totsum(minind)/dofsub,'.',...
       'markersize',20,'color',get_color(ifield),'HandleVisibility','off');
chi2dat.bestg1(ifield).g1=G1_arr(minind);
end
plot(G1_arr,subchi2_allfield/subdof_allfield,'.-','color',...
    get_color(ifield+5),'markersize',20,'DisplayName','all fields');hold on
minind=find(subchi2_allfield==min(subchi2_allfield));
plot(G1_arr(minind),subchi2_allfield(minind)/subdof_allfield,'*',...
    'markersize',20,'color',get_color(ifield+5),'HandleVisibility','off');

if inst==1
    ylim([0,50]);
else
    ylim([0,300]);
end
yL = get(gca,'YLim');
line([G1_FF G1_FF],yL,'linestyle','--','Color',get_color(3),...
    'DisplayName','FFlab');
legend('show','location','northwest')
legend boxoff
xlabel('$G1$','interpreter','latex','fontsize',18)
ylabel('$\chi^2/dof$',...
    'interpreter','latex','fontsize',18)

savename=strcat(savedir,'chi2g1_fit');
print(savename,'-dpng');%close

minind=find(subchi2_allfield==min(subchi2_allfield));
chi2dat.bestg1_allfields=G1_arr(minind);
save(strcat(savedir,'chi2dat'),'chi2dat');
return