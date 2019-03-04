function r_arr = get_mask_radius(inst,ifield,m_arr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Masking function
% Input:
% m_arr - I/H band magnitude from catalog.
% r_arr - masking radius in arcsec
%
% r_arr is set to 3.5 if it is <3.5 from the fitting funciton
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if (inst ==1) && (ifield == 4)
    p = [4.35614072e-04 -3.42183125e-02 1.08173147e+00 -1.64779272e+01 1.02339117e+02];
elseif (inst ==1) && (ifield == 5)
    p = [3.52991997e-03 -2.49966754e-01 6.69611761e+00 -8.12775345e+01 3.82246341e+02];
    % test aggresive masking
    %p = [5.23791526e-05 -7.57545501e-02 4.27210967e+00 -8.06530129e+01 5.05553967e+02];
elseif (inst ==1) && (ifield == 6)
    p = [2.73314825e-03 -2.02562607e-01 5.68442521e+00 -7.21790791e+01 3.53582782e+02];
elseif (inst ==1) && (ifield == 7)
    p = [1.42515137e-03 -1.16957156e-01 3.65605478e+00 -5.16458485e+01 2.79057652e+02];
elseif (inst ==1) && (ifield == 8)
    p = [3.37845041e-04 -2.78816481e-02 9.01832457e-01 -1.37387245e+01 8.42207024e+01];
elseif (inst ==2) && (ifield == 4)
    p = [5.53301366e-03 -3.37451258e-01 7.94797962e+00 -8.73119102e+01 3.84623370e+02];
elseif (inst ==2) && (ifield == 5)
    p = [3.62033577e-03 -2.49849600e-01 6.56485125e+00 -7.87194708e+01 3.68374414e+02];
elseif (inst ==2) && (ifield == 6)
    p = [4.92069543e-03 -3.18837417e-01 7.91316610e+00 -9.03031944e+01 4.05781437e+02];
elseif (inst ==2) && (ifield == 7)
    p = [1.99557798e-03 -1.47949923e-01 4.26889685e+00 -5.68580896e+01 2.95139426e+02];
elseif (inst ==2) && (ifield == 8)
    p = [2.98636348e-03 -1.98612472e-01 5.12166410e+00 -6.13109340e+01 2.90143425e+02];

end

% set the radius to 0.5 pix if the function goes below 0.5
mbins = min(m_arr):0.1:max(m_arr);
r_arr = polyval(p,mbins);
mcut = find(r_arr < 0.5);
if numel(mcut) > 0
    mcut = mbins(mcut(1));
else
    mcut = max(mbins) + 1;
end

% set the radius to min point if the function starts increase with mag
r_arr = polyval(p,m_arr);
r_arr(m_arr > mcut) = 0.5;

sp = find(r_arr == min(r_arr));
if numel(sp)>=2
    sp = sp(1);
end
r_arr(m_arr > m_arr(sp)) = r_arr(sp);

r_arr = r_arr.*7;



%     if ifield==4
%         a1 = 754.2;
%         b1 = -3.637;
%         c1 = 11.11;
%         r_arr = a1.*exp(-((m_arr-b1)./c1).^2);
%     elseif ifield==5
%         a1 = 1.417e+04;
%         b1 = -25.49;
%         c1 = 16.87;
%         r_arr = a1.*exp(-((m_arr-b1)./c1).^2);
%     elseif ifield==6
%         a0 = 180.1;
%         a1 = 218.2;
%         b1 = 85.37;
%         a2 = 51.56;
%         b2 = 36.7;
%         w = 0.1714;
%         r_arr = a0 + a1.*cos(m_arr.*w) + b1.*sin(m_arr.*w) + ...
%                a2.*cos(2.*m_arr.*w) + b2.*sin(2.*m_arr.*w);
%     elseif ifield==7
%         p1 = -0.2016;
%         p2 = 11.55;
%         p3 = -222.9;
%         p4 = 1458;
%         r_arr = p1.*m_arr.^3 + p2.*m_arr.^2 + p3.*m_arr + p4;
%     elseif ifield==8
%         r_arr = 208.4.*exp(-(m_arr-3.082).^2./(9.209).^2) + ...
%             31.83.*exp(-(m_arr-9.297).^2./(1.731).^2); 
%     end
% 
% r_arr(r_arr<3.5) = 3.5;

return