function [m_arr, I_arr] = get_corrected_mag(inst, mlin_arr, my_arr, cls_arr)
% Get the corrected I, H band magnitude for PanSTARRS catalog
% The correction factors are the results from mag_correction

sr = ((7./3600.0)*(pi/180.0)).^2;
    if inst==1
        lambdaeff=1.05;
        Ilin_arr=3631*10.^(-mlin_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
        I_arr = 1.07 * Ilin_arr;
        sp = find(my_arr < 23.3 & cls_arr==1);
        I_arr(sp) = (-0.44 * my_arr(sp) + 10.34) .* Ilin_arr(sp);
        sp = find(my_arr >= 23.3 & cls_arr==1);
        I_arr(sp) = (-0.44 * 23.3 + 10.34) .* Ilin_arr(sp);
    else
        lambdaeff=1.79;
        Ilin_arr=3631*10.^(-mlin_arr/2.5)*(3/lambdaeff)*1e6/(sr*1e9);
        I_arr = 1.40 * Ilin_arr;
        sp = find(my_arr < 23.3 & cls_arr==1);
        I_arr(sp) = (-0.47 * my_arr(sp) + 11.11) .* Ilin_arr(sp);
        sp = find(my_arr >= 23.3 & cls_arr==1);
        I_arr(sp) = (-0.47 * 23.3 + 11.11) .* Ilin_arr(sp);
    end
    
    m_arr = mlin_arr - 2.5 * log10(I_arr./Ilin_arr);

return