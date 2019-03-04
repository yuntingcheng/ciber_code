function [map,offset] = fastlinefit_frin(frames,quadrant,nskip,nframes)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% works like reduc_fit_wrapper only much dumber and simpler.
% should only be used to make simple maps, runs about 25% faster with no
% bells or whistles.
% This is mostly useful for things like focus data which 
% need to be turned around quick
%
% To run this you need to specify a quadrant.  If you want the whole map it
% needs a 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% stuff for troubleshooting 
%fitparams = get_default_fitparams;
%fitparams.labdata = 1;
%field = '12-12-2011_10-19-21';
%inst =4;
%fitparams.nskip =8;
 %% Parse input arguments
 
%% set up map parameters
s = size(frames);
map = zeros(s(2));
x = 1:s(1);
offset = zeros(s(2));

%% Set up you bounds based on what quadrant you're in

if quadrant == 0
    botc = 1;
    topc = s(2);
    botr = 1;
    topr = s(3);
end

if quadrant == 1
    botc = 1;
    topc = s(2)/2;
    botr = 1;
    topr = s(3)/2;
end

if quadrant == 2
    botc = 1;
    topc = s(2)/2;
    botr = s(3)/2;
    topr = s(3);
end

if quadrant == 3
    botc = s(2)/2;
    topc = s(2);
    botr = s(3)/2;
    topr = s(3);
end

if quadrant == 4
    botc = s(2)/2;
    topc = s(2);
    botr = 1;
    topr = s(2)/2;
end

%% loop around array and fit lines

%for cc=botc:topc
%    for rr=botr:topr       
%        line = linfit(x(:),squeeze(frames(:,cc,rr)));       
%        map(cc,rr) = line(1);
%    end
%    if cc/20 == round(cc/20) || cc ==topc
%        display(strcat(num2str(round(100.*(cc - botc)./(topc - botc))),'% done with line fit') )
%    end
%end

yp = squeeze(frames(:,20,20));
%figure(1)
%plot(yp)


for cc=botc:topc
    for rr=botr:topr       
        xp = x(:);
        yp = squeeze(frames(:,cc,rr));
        Vp(:,2) = ones(length(xp),1,class(xp));
        Vp(:,1) = xp.*Vp(:,2);
        [Qp,Rp] = qr(Vp,0);
        pp = Rp\(Qp'*yp);    % Same as p = V\y;
        line = pp.';
        %line = linfit(x(:),squeeze(frames(:,cc,rr)));       
        map(cc,rr) = line(1);
        offset(cc,rr) = line(2);
    end
    if cc/20 == round(cc/20) || cc ==topc
        %display(strcat(num2str(round(100.*(cc - botc)./(topc - botc))),'% done with line fit') )
    end
end

%if inst ==1 || inst ==2
%    map = map.*(-1);
%end

%figure(2)
%imageclip(map);
%set(gca,'FontSize',16)
%title(strcat(field,'Dumb Simple Map'))

return