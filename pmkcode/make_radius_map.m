function rad = make_radius_map(map,cenx,ceny)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Given an input map and a specified center, this will
% create a map with each pixels value its distance from
% the specified pixel.
% 
% example
% 
%
% rad = make_radius_map(map,120,120)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


rad = map.*0;
s=size(rad);
for i=1:s(1)
	for ii=1:s(2)
		 rad(i,ii)=sqrt( (cenx-i)^2 + (ceny -ii)^2  );
    end
end

return
