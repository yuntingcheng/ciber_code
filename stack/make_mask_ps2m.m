function mask = make_mask_ps2m(flight,inst,ifield,m_min,m_max,m_join)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Produce the mask from PanSTARR, and join with 2MASS mask below m_join.
% (PanSTARRS not complete at bright end)
%Input:
%(Reqiured)
% - flight: flight # (40030 for 4th flight)
% - inst: 1 or 2 (I/H)
% - ifield: 4,5,6,7,8 
% - m_min: min masking magnitude (PS y band)
% - m_max: max masking magnitude (PS y band)
% - m_join: mag below which use both catalog
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mask = make_mask_ps(flight,inst,ifield,0,m_min,m_max);

if m_min<m_join
    tmmask = make_mask_2m(flight,inst,ifield,m_min,m_join);
    mask = mask.*tmmask;
end

end
