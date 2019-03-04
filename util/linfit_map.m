function [slopedata offsetdata vardata]= linfit_map(frames,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do line fit to raw time-stream data with linfit (MZ code)
% Copy from see reduc_line_fit line 260
% Remove the NaN before fitting
%
% fast=1: use method in fastlinefit_frin.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Parse data
  p = inputParser;
  p.addRequired('frames',@isnumeric);
  p.addOptional('verbose',0,@isnumeric);
  p.addOptional('fast',1,@isnumeric);
  
  p.parse(frames,varargin{:});
  frames     = p.Results.frames;
  verbose  = p.Results.verbose;
  fast = p.Results.fast;
  
  clear p varargin;
%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,Nx,Ny]=size(frames);
slopedata=zeros(size(frames,2),size(frames,3));
offsetdata=zeros(size(frames,2),size(frames,3));
vardata=zeros(size(frames,2),size(frames,3));

for ind=1:Nx
    for jnd=1:Ny
        yp=frames(:,ind,jnd);
        xp=1:numel(yp);

        %remove the NaN in fitting
        xp(find(yp~=yp))=[];
        yp(find(yp~=yp))=[];

        if fast==0
            % do linfit
            [thisfit,fitcov] = linfit(xp,yp);
            slopedata(ind,jnd)  = thisfit(1);
            offsetdata(ind,jnd) = thisfit(2);
            vardata(ind,jnd)    = fitcov(1,1);
        end

        if fast~=0
            Vp(:,2) = ones(length(xp),1,class(xp));
            Vp(:,1) = xp(:).*Vp(:,2);
            [Qp,Rp] = qr(Vp,0);
            pp = Rp\(Qp'*yp);    % Same as p = V\y;
            line = pp.';
            slopedata(ind,jnd) = line(1);
            offsetdata(ind,jnd) = line(2);
            vardata(ind,jnd) = 0; % no covariance return!!
        end

    end

    % print out the progress ~every 10%
    if verbose
        if mod(ind,round(Nx/10))==0
        pr=sprintf('linfit_map done %d/%d=%.1f%%',...
            ind,Nx,ind/Nx*100);disp(pr);
        end
    end
end
return