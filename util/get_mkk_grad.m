function [Mkk] = get_mkk_grad(mask,pixsize,lbin,numsims,logbins,w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get_mkk_grad.m - Get the mode-coupling correction matrix w/ simulations
%
% For each simulation, we remove a gradient w/ plane_fit.m
%
%  Inputs:
%    - mask:  The mask for this particular Mkk.
%    - pixsize: The pixel scale of the mask in arcsec.
%    - lbin: The l-mode bin edges.  
%    - numsims: Number of simulations to calculate Mkk.  
%    - logbins(0/1)
%    - w: weight in 2D Fourier space
%
% Outputs:
%   - Mkk:  The mode-coupling correction matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = numel(lbin) - 1;
Mkk = zeros(M);
[dim1 dim2] = size(mask);

% Build array to hold N power spectra used to calculate Mkk.
Cll = zeros(numsims,M);

for j = 1:M
    Cl1 = zeros(1,M);
      for i=1:numsims
        disp(sprintf('At ell = %d, numsims=%d\n',j,i));
        Cl1(j) = 1;
        map = map_from_power_spec(lbin,Cl1,dim1,dim2,pixsize,1);
        [Cl0] = get_angular_spec(map,map,pixsize,'logbins',logbins);
        Clnorm=1;
        if Cl0(j)~=0
            Clnorm=Cl0(j);
        end
        
        if Cl0(j)==0
            fit_map=zeros(size(mask));
        else
            [fit_map]=plane_fit(map,mask);
            fit_map=fit_map-mean(fit_map(find(fit_map)));
        end
        
        mfm=(map-fit_map).*mask;mfm=mfm-mean(mfm(find(mfm)));mfm=mfm.*mask;

        % Calculate the power spectrum which becomes the ith row of Mkk
        [~,~,~,~,~,~,Cl2d] = get_angular_spec(mfm,mfm,pixsize);
        Cl=Cl_from_Cl2d(Cl2d,pixsize,'w',w);
        
        if Cl0(j)~=0
            Clnorm=Cl0(j);Cl=Cl./Clnorm;
        end
        
        Cll(i,:) = Cl;
      end
    % The average of M simulations is the resulting row.
    Mkk(j,:) = squeeze(nanmean(Cll));
end
Mkk = Mkk';
end