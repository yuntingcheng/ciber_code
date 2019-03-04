dd=('/Users/ytcheng/ciber/data/40030/framedata/');

files = dir(strcat(dd,'TM1*'));

close = (1:numel(files))*0;
t1 = close;
time = t1;
op = close;
lamp = close;
ts=close;
for i=1:numel(files)
    i/numel(files)
    load(strcat(dd,files(i).name));
    t1(i) = HK.T1;
    time(i) = HK.ttt;
    close(i) = HK.close;
    op(i) = HK.open;
    lamp(i) = HK.lamp;
    ts(i) = arraymap(821,168);

end
%%
%mint = (time-time(1))/60 - 3.542;
plot(mint,close)
hold on
plot(mint,t1,'color','red')
plot(mint,op,'color','green')
hold off
%axis([-6,20,-1,2])
%%
l=-4;
h=5;

subplot(3,1,1)
plot(mint,ts);
axis([l,h,3250,3350])

subplot(3,1,2)
plot(mint,t1)
axis([l,h,-1,2])

subplot(3,1,3)
plot(mint,op,'color','red');
hold on
plot(mint,close,'color','black');
plot(mint,lamp,'color','green');

hold off
axis([l,h,-1,2])
%%

nfr=numel(files);
label_arr=1:nfr;
subplot(3,1,1)
plot(label_arr,ts);
axis([0,600,3250,3350])

subplot(3,1,2)
plot(label_arr,t1)
axis([0,600,-1,2])

subplot(3,1,3)
plot(label_arr,op,'color','red');
hold on
plot(label_arr,close,'color','black');
plot(label_arr,lamp,'color','green');

hold off
axis([0,600,-1,2])
%%
figure
tl=97;
th=119;
nfr=numel(files);
label_arr=1:nfr;
subplot(3,1,1)
plot(label_arr(tl:th),ts(tl:th));
axis([tl-10,th+10,3250,3350])

subplot(3,1,2)
plot(label_arr(tl:th),t1(tl:th))
axis([tl-10,th+10,-1,2])

subplot(3,1,3)
plot(label_arr(tl:th),op(tl:th),'color','red');
hold on
plot(label_arr(tl:th),close(tl:th),'color','black');
plot(label_arr(tl:th),lamp(tl:th),'color','green');

hold off
axis([tl-10,th+10,-1,2])
