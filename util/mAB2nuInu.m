function I=mAB2nuInu(m,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Convert AB mag to CIBER surface brightness
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% CIBER effective wavelength 
Ilambda=1.05; % [um] 
Hlambda=1.79; % [um] 
nuI=3e8/(Ilambda*1e-6); % [Hz]
nuH=3e8/(Hlambda*1e-6); % [Hz]

if inst==1
    nu=nuI;
else
    nu=nuH;
end

%%% Omega_pix %%%%%
sr=(7./180/60/60*pi)^2; % [Sr]

%%% convert unit
F=3631.*10.^(m/-2.5); %[Jy]
I=F*10^-17/sr*nu;

return
