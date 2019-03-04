function cal = get_cal_apf2nWpm2ps(inst)
% cal factor from ADU/fr to nW/m2/sr. 
% The values are the results from fit_calfac.m

if inst ==1
    cal(4).apf2nWpm2ps = -363.60;
    cal(5).apf2nWpm2ps = -308.59;
    cal(6).apf2nWpm2ps = -349.96;
    cal(7).apf2nWpm2ps = -336.06;
    cal(8).apf2nWpm2ps = -310.44;
else
    cal(4).apf2nWpm2ps = -119.32;
    cal(5).apf2nWpm2ps = -115.59;
    cal(6).apf2nWpm2ps = -119.50;
    cal(7).apf2nWpm2ps = -121.67;
    cal(8).apf2nWpm2ps = -117.65;    
end

return