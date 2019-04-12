function cal = get_cal_apf2nWpm2ps(inst)
% cal factor from ADU/fr to nW/m2/sr. 
% The values are the results from fit_calfac.m

if inst ==1
    cal(4).apf2nWpm2ps = -347.92;
    cal(5).apf2nWpm2ps = -305.81;
    cal(6).apf2nWpm2ps = -369.32;
    cal(7).apf2nWpm2ps = -333.67;
    cal(8).apf2nWpm2ps = -314.33;
else
    cal(4).apf2nWpm2ps = -117.69;
    cal(5).apf2nWpm2ps = -116.20;
    cal(6).apf2nWpm2ps = -118.79;
    cal(7).apf2nWpm2ps = -127.43;
    cal(8).apf2nWpm2ps = -117.96;    
end

return