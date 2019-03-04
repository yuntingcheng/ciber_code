function [shotmap]=photonnoise_realization(slopemap,g1,nfr,frate)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function generate photon noise realization given 
%raw slope map and gain factor g1,g2.
%Ref: get_shotmap.m
%
%Input:
%   -slopemap: raw slope map(ADU/fr)
%   -g1:gain factor g1 [(e-/s)/(ADU/fr)]
%   -nfr:# of frames in slopemap
%   -frate:frame rate
%Output:
%   -shotmap:sim photon noise map (ADU/fr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % convert to eps
    flightsignal = slopemap.*g1;
    
    %Make shot-noise term from flight data Garnett and Forrest
    shot_sigma = sqrt((6./5).*(flightsignal./(nfr./frate))...
            .*((nfr.^2 +1)/(nfr.^2 -1)));

    shotmapeps=normrnd(zeros(size(slopemap)),real(shot_sigma));
    
    shotmap=shotmapeps./g1;
end
