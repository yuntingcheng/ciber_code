function [stampercb,stamperps,maskstamper]=...
    stackihl_ps_randomN(N_arr,dx,cbmap,psmap,mask,savename,iter)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Stack src based on PanSTARRS catalog
%
%Input:
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - type: 1:gal, -1:star, 0:all, 2:undefined
% - m_min: min masking magnitude
% - m_max: max masking magnitude
% - dx: stamp size is 2*dx+1, dx is 10 times finer CIBER pix
% - cbmap:
% - psmap:
% - mask_inst:
% - strmask:
% - strnum:
% - verbose:1 or 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% randomize the position
subx_arr=1023*rand(1,max(N_arr))+1;
suby_arr=1023*rand(1,max(N_arr))+1;

stampercb=zeros(2*dx+1);
stamperps=zeros(2*dx+1);
maskstamper=zeros(2*dx+1);

for i=1:max(N_arr)
    cbmapi = cbmap.*mask;
    psmapi = psmap.*mask;
    %%% zero padding
    cbmapi = padarray(cbmapi,[dx/10 dx/10],0,'both');
    psmapi = padarray(psmapi,[dx/10 dx/10],0,'both');
    %%% rebin to 10x finer map
    stampcb = imresize(cbmapi,10,'nearest');
    stampps = imresize(psmapi,10,'nearest');
    %%% get the stamp
    xcent = round(subx_arr(i)*10-4.5) + dx;
    ycent = round(suby_arr(i)*10-4.5) + dx;
    stampcb = stampcb(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    stampps = stampps(xcent-dx:xcent+dx,ycent-dx:ycent+dx);
    maskstamp=zeros(size(stampcb));
    maskstamp(find(stampcb))=1;
    %%% stack
    if mod(i,4)==0
        stampercb=stampercb+imrotate(stampcb, 0);
        stamperps=stamperps+imrotate(stampps, 0);
        maskstamper=maskstamper+imrotate(maskstamp, 0);
    elseif mod(i,4)==1
        stampercb=stampercb+imrotate(stampcb, 90);
        stamperps=stamperps+imrotate(stampps, 90);
        maskstamper=maskstamper+imrotate(maskstamp, 90);
    elseif mod(i,4)==2
        stampercb=stampercb+imrotate(stampcb, 180);
        stamperps=stamperps+imrotate(stampps, 180);
        maskstamper=maskstamper+imrotate(maskstamp, 180);
    elseif mod(i,4)==3
        stampercb=stampercb+imrotate(stampcb, 270);
        stamperps=stamperps+imrotate(stampps, 270);
        maskstamper=maskstamper+imrotate(maskstamp, 270);
    end

    if ismember(i,N_arr)
%         profile = radial_prof(stampercb./maskstamper,ones(2*dx+1),...
%             dx+1,dx+1,1,25,'sig',3,'iter_clip',3);
        profile = radial_prof(stampercb./maskstamper,ones(2*dx+1),...
            dx+1,dx+1,1,25,'weight',maskstamper);
        r_arr=profile.r*0.7;
        profcb_arr=(profile.prof);
        errcb_arr=profile.err;
        
%         profile = radial_prof(stamperps./maskstamper,ones(2*dx+1),...
%             dx+1,dx+1,1,25,'sig',3,'iter_clip',3);
        profile = radial_prof(stamperps./maskstamper,ones(2*dx+1),...
            dx+1,dx+1,1,25,'weight',maskstamper);
        profps_arr=(profile.prof);
        errps_arr=profile.err;

        M(:,1) = r_arr;
        M(:,2) = profcb_arr;
        M(:,3) = errcb_arr;
        M(:,4) = profps_arr;
        M(:,5) = errps_arr;
        csvwrite(strcat(savename,'_randprof',...
            '_',num2str(i),'_',num2str(iter),'.csv'),M);
        
        %{
        fits_write(strcat(savename,'_stampercb',...
            '_',num2str(i),'_',num2str(iter),'.fits'),stampercb);
        fits_write(strcat(savename,'_stamperps',...
            '_',num2str(i),'_',num2str(iter),'.fits'),stamperps);
        fits_write(strcat(savename,'_maskstamper',...
            '_',num2str(i),'_',num2str(iter),'.fits'),maskstamper);
        %}
        
        disp(sprintf('stack %d random src',i));
    end
      
end

end
