function c=get_color(ind)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is a list of distinguishable color copy from
% http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/
% The function return the color code of i-th color in the list
% example:
% plot(x_arr,y_arr,'color',get_color(ind))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colororder = [
	0.00  0.00  1.00
	0.00  0.50  0.00 
	1.00  0.00  0.00 
	0.00  0.75  0.75
	0.75  0.00  0.75
	0.75  0.75  0.00 
	0.25  0.25  0.25
	0.75  0.25  0.25
	0.95  0.95  0.00 
	0.25  0.25  0.75
	0.75  0.75  0.75
	0.00  1.00  0.00 
	0.76  0.57  0.17
	0.54  0.63  0.22
	0.34  0.57  0.92
	1.00  0.10  0.60
	0.88  0.75  0.73
	0.10  0.49  0.47
	0.66  0.34  0.65
	0.99  0.41  0.23
];
c=colororder(rem(ind,20)+1,:);
return
