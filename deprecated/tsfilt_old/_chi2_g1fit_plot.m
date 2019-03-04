%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plot the Chi2-G1 of each field, and chi2 fit with 
%all the fields.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flight=40030;
inst=1;
fields_info=get_fields(flight,inst);
savedir='/Users/ytcheng/ciber/doc/20170209_TsFilter/G1fit/';
%load(strcat(savedir,'chi2dat'),'chi2dat');
load(strcat(savedir,'chi2datw'),'chi2dat');
G1_arr=-4.5:0.1:-2.5;
nfr_arr=2:19;
figure
subchi2allfield=zeros(size(G1_arr));
subdofallfield=-1;
bestg1=zeros(1,9);
for ifield=[1,2,4,5,6,7,8]
name=fields_info(ifield).name;
%figure
chi2totsum=zeros(1,numel(G1_arr));
subchi2totsum=zeros(1,numel(G1_arr));
doftot=-1;
dofsub=-1;
for nfr=nfr_arr
    field_use=chi2dat(nfr).field_use;
    chi2tot_arr=chi2dat(nfr).chi2tot_arr;
    subchi2tot_arr=chi2dat(nfr).subchi2tot_arr;
    if sum(find(field_use==ifield))
        ind=find(field_use==ifield);
        chi2tot=chi2tot_arr(ind,:);
        subchi2tot=subchi2tot_arr(ind,:);
        chi2totsum=chi2totsum+chi2tot;
        subchi2totsum=subchi2totsum+subchi2tot;
        doftot=doftot+21;
        dofsub=dofsub+5;
        %plot(G1_arr,chi2tot,'o--','color',get_color(nfr-1));hold on
        %plot(G1_arr,subchi2tot,'o-','color',get_color(nfr-1));hold on
    end  
end
if find(ifield==[4,5,6,7,8])
    subchi2allfield=subchi2allfield+subchi2totsum;
    subdofallfield=subdofallfield+dofsub;
end
%figure
%plot(G1_arr,chi2totsum/doftot,'o-');hold on
plotdat(ifield).subchi2totsum=subchi2totsum;
plot(G1_arr,subchi2totsum/dofsub,'o-','color',get_color(ifield),...
    'DisplayName',name);hold on
minind=find(subchi2totsum==min(subchi2totsum));
plot(G1_arr(minind),subchi2totsum(minind)/dofsub,'.',...
       'markersize',20,'color',get_color(ifield),'HandleVisibility','off');
%title(sprintf('field%d',ifield));
bestg1(ifield)=G1_arr(minind);
end
yL = get(gca,'YLim');
plotdat(9).subchi2totsum=subchi2totsum;
plot(G1_arr,subchi2allfield/subdofallfield,'.-','color',...
    get_color(ifield+5),'markersize',20,'DisplayName','all fields');hold on
line([-2.68 -2.68],yL,'linestyle','--','Color',get_color(3),...
    'DisplayName','FFlab');
legend('show','location','southeast')
legend boxoff
ylim([0,20])
xlabel('$G1$','interpreter','latex','fontsize',18)
ylabel('$\chi^2/dof$',...
    'interpreter','latex','fontsize',18)

%savename=strcat(savedir,'chi2g1');
savename=strcat(savedir,'chi2gw');
print(savename,'-dpng');%close

minind=find(subchi2allfield==min(subchi2allfield));
bestg1(9)=G1_arr(minind);
%save(strcat(savedir,'bestg1'),'bestg1');
save(strcat(savedir,'bestg1w'),'bestg1');