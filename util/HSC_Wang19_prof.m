function [I_arr, I_arr_sub,R200,param_arr] = HSC_Wang19_prof(r_arr,im,extendedness)
% https://arxiv.org/pdf/1811.04714.pdf table 3
% r_arr [arcsec]
% im=1:4 - 16-17;17-18;18-19;19-20
% extendedness: true - C > 2.6; false - C < 2.6

% median R200 from SIDES abundnace matching
R200_arr = [98.90,62.83,42.48,29.34]; % [arcsec]
R200 = R200_arr(im);

if extendedness
    param1 = [-8.471,1.5320,0.0056];
    param2 = [-8.9330,2.6190,0.0165];
    I1_arr = sersic(r_arr./R200,param1(1),param1(2),param1(3));
    I2_arr = sersic(r_arr./R200,param2(1),param2(2),param1(3)+param2(3));
    I_arr = I1_arr + I2_arr;
    I_arr_sub = zeros([2,size(I_arr)]);
    I_arr_sub(1,:,:) = I1_arr;
    I_arr_sub(2,:,:) = I2_arr;
    param_arr = [param1;param2];
    
else
    param1 = [-7.3500,0.0101,0.0015];
    param2 = [-8.7930,1.143,0.0231];
    param3 = [-9.9220,3.691,0.0001];
    I1_arr = sersic(r_arr./R200,param1(1),param1(2),param1(3));% This 0.0101 is weird
    I2_arr = sersic(r_arr./R200,param2(1),param2(2),param1(3)+param2(3));
    I3_arr = sersic(r_arr./R200,param2(1),param2(2),param1(3)+param2(3)+param3(3));
    I_arr =  I1_arr + I2_arr + I3_arr;
    I_arr_sub = zeros([3,size(I_arr)]);
    I_arr_sub(1,:,:) = I1_arr;
    I_arr_sub(2,:,:) = I2_arr;
    I_arr_sub(3,:,:) = I3_arr;
    param_arr = [param1;param2;param3];

end

return