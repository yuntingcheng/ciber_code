function [name,mask] = HSC_fields_info(idx)
%%%%%%%%%%%%%%
% recording the processed HSC field name and their index
%%%%%%%%%%%%%%
%%
mask = ones(1024);
name = strcat('W05_',num2str(idx));

if idx==0
    mask=circular_mask(622,909,20,mask);
    mask=circular_mask(279,915,30,mask);
    mask=circular_mask(527,307,25,mask);
    mask=circular_mask(380,788,20,mask);
    mask=circular_mask(430,191,15,mask);
elseif idx==1
    mask=circular_mask(633,252,20,mask);
    mask=circular_mask(812,225,20,mask);
elseif idx==2
    mask=circular_mask(641,879,20,mask);
    mask=circular_mask(134,369,20,mask);
    mask=circular_mask(222,679,20,mask);
elseif idx==3
    mask=circular_mask(629,384,45,mask);
    mask=circular_mask(790,219,25,mask);
elseif idx==4
    mask=circular_mask(242,710,60,mask);
    mask=circular_mask(963,266,40,mask);
elseif idx==5
    mask=circular_mask(294,508,20,mask);
elseif idx==6
    mask=circular_mask(196,535,60,mask);
    mask=circular_mask(95,290,25,mask);
    mask=circular_mask(562,435,50,mask);
    mask=circular_mask(870,314,45,mask);
    mask=circular_mask(823,189,20,mask);
elseif idx==7
    mask=circular_mask(183,331,15,mask);
    mask=circular_mask(217,779,20,mask);
    mask=circular_mask(554,150,45,mask);
    mask=circular_mask(638,351,25,mask);
    mask=circular_mask(774,191,50,mask);
    mask=circular_mask(905,100,30,mask);
    mask=circular_mask(1008,938,50,mask);
    mask=circular_mask(946,934,10,mask);
    mask=circular_mask(508,419,25,mask);
    mask=circular_mask(267,97,15,mask);
elseif idx==8
    mask=circular_mask(291,452,20,mask);
    mask=circular_mask(462,483,20,mask);
    mask=circular_mask(437,508,15,mask);
    mask=circular_mask(880,484,25,mask);
elseif idx==9
    mask=circular_mask(128,792,25,mask);
    mask=circular_mask(239,973,25,mask);
    mask=circular_mask(894,77,25,mask);
elseif idx==10
    mask=circular_mask(185,346,25,mask);
    mask=circular_mask(467,905,60,mask);
    mask=circular_mask(625,970,40,mask);
    mask=circular_mask(815,990,25,mask);
elseif idx==11
    mask=circular_mask(175,856,50,mask);
    mask=circular_mask(701,517,20,mask);
    mask=circular_mask(914,710,15,mask);
    mask=circular_mask(150,21,15,mask);


end

%%
% mask = ones(642);
% if idx==1
%     name = 'UD_COSMOS';
% elseif idx==2
%     name = 'D_DEEP2-3';
% elseif idx==3
%     name = 'D_ELAISN1';
% 
% elseif idx==4
%     name = 'XMM_00';
% elseif idx==5
%     name = 'XMM_01';
% elseif idx==6
%     name = 'XMM_10';
% elseif idx==7
%     name = 'XMM_11';
% elseif idx==8
%     name = 'XMM_20';
% elseif idx==9
%     name = 'XMM_21';
% elseif idx==10
%     name = 'XMM_30';
% elseif idx==11
%     name = 'XMM_31';
% 
% elseif idx==12
%     name = 'GAMA09H_00';
% elseif idx==13
%     name = 'GAMA09H_01';
% elseif idx==14
%     name = 'GAMA09H_10';
% elseif idx==15
%     name = 'GAMA09H_11';
% elseif idx==16
%     name = 'GAMA09H_20';
% elseif idx==17
%     name = 'GAMA09H_21';
% elseif idx==18
%     name = 'GAMA09H_30';
% elseif idx==19
%     name = 'GAMA09H_31';
% elseif idx==20
%     name = 'GAMA09H_40';
% elseif idx==21
%     name = 'GAMA09H_41';
% 
% end
return