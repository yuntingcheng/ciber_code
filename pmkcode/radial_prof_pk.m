function profile = radial_prof_pk(map,mask,cenx,ceny,nbins)

  %% make a map containing the distance
  % of each pixel to the desired center

  rad = make_radius_map(map,cenx,ceny);

  %% setup output parameters

  r=1:nbins;
  prof=r;
  err=r;
  prof2d = rad.*0;
  prof2d_unholy = rad.*0;
  
  npix=r;
  %% Get the limits for the loop

  topw = max(rad(:));
  botw = 0;
  dw = (topw - botw)/nbins;

  %% loop around the distance map and bin up your signal map

  for i=1:nbins
      spot=find(rad >= ((i-1)*dw + botw) & rad < ((i)*dw + botw) & mask);
      spotnom=find(rad >= ((i-1)*dw + botw) & rad < ((i)*dw + botw));
      vecr=rad(spot); 
      r(i)=nanmean(vecr(:));
      vecim=map(spot);
      prof(i)=nanmean(vecim(:));
      prof2d(spot) = prof(i);
      prof2d_unholy(spotnom) = prof(i);
      %err(i)=std(vecim(:))  /   sqrt(  length(vecim(:))  );
      err(i)=nanmean(vecim(:))  /   sqrt(  length(vecim(:))  );
      npix(i)=length(vecim(:));
   end

   %% make plots
   %{
   fig=figure;
   plot(r,prof,':')
   errorbar(r,prof,err)
   ax = get(fig,'CurrentAxes');
   set(ax,'XScale','log','YScale','log');
   title('Radial Profile','FontSize',16)
   xlabel('Radial Distance (Pixels)','FontSize',16)
   ylabel('Amplitude','FontSize',16)
   %}
   %% get a nice structure to pass through to matlab
  
   profile = struct('r',r,'prof',prof,'err',err,...
       'prof2d',prof2d,'prof2d_unholy',prof2d_unholy);


return

