function wl = bands_wavelengths(name,band)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The band wavelength of surveys in um.
% name: 'CIBER', '2MASS', 'UKIDSS', 'PanSTARRS'
% band: 
% - CIBER: 'I', 'H'
% - 2MASS: 'j', 'h', 'k'
% - UKIDSS: 'y', 'j', 'h', 'k'
% - PanSTARRS: 'g', 'r', 'i', 'z', 'y'
%
% example:
% wl = band_wavelengths('PanSTARRS', 'y')
% wl
% >> 0.9633
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(name,'CIBER')
    keySet =   {'I','H'};
    if ~ismember(band,keySet)
        disp(sprintf('survey name not exist!!!'));
        wl=nan;
    else
    wlSet = [1.05, 1.79];
    wlObj = containers.Map(keySet,wlSet);
    wl = wlObj(band);
    end
elseif strcmp(name,'2MASS')
    keySet =   {'j','h','k'};
    if ~ismember(band,keySet)
        disp(sprintf('survey name not exist!!!'));
        wl=nan;
    else
    wlObj = containers.Map(keySet,wlSet);
    wl = wlObj(band);
    end
elseif strcmp(name,'UKIDSS')
    keySet =   {'y','j','h','k'};
    if ~ismember(band,keySet)
        disp(sprintf('survey name not exist!!!'));
        wl=nan;
    else
    wlSet = [1.0305,1.2483,1.6313,2.2010];
    wlObj = containers.Map(keySet,wlSet);
    wl = wlObj(band);
    end
elseif strcmp(name,'PanSTARRS')
    keySet =   {'g','r','i','z','y'};
    if ~ismember(band,keySet)
        disp(sprintf('survey name not exist!!!'));
        wl=nan;
    else
    wlSet = [0.4866,0.6215,0.7545,0.8679,0.9633];
    wlObj = containers.Map(keySet,wlSet);
    wl = wlObj(band);
    end
else
    disp(sprintf('survey name not exist!!!'));
    wl=nan;
end    
    
    

end
