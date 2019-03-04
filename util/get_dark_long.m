function [time_arr,start_arr,end_arr]=get_dark_long(flight,band)

if flight==36265

    switch band
    case 1
    time_arr={'12-35-09';'18-10-47';'18-18-01';'18-18-01';'20-00-01';...
              '21-11-35'};
    start_arr=[3,3,20,85,22,3];
    end_arr=[39,65,82,129,48,62];
    case 2
    time_arr={'12-35-09';'18-10-47';'18-18-01';'18-18-01';'20-00-01';...
              '21-11-35'};
    start_arr=[3,3,19,84,3,3];
    end_arr=[39,65,81,129,48,63];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif flight==36277
    switch band
    case 1
    time_arr={'12-12-17';'12-55-11';'12-55-11';'13-01-18';...
              '13-01-18';'13-05-21';'13-05-21';'13-10-11';'13-10-11';...
              '13-14-28';'13-17-00';'13-17-00';'08-24-00';'08-29-36';...
              '08-58-51';'09-51-30';'10-07-23';'10-32-16';'10-58-11';...
              '14-40-01';'14-49-10'};
    start_arr=[3,3,65,3,54,3,48,16,81,3,3,64,3,8,3,3,3,3,3,6,23];
    end_arr=[50,62,127,51,88,45,110,78,110,57,61,118,...
                                               36,35,34,38,43,54,29,63,66];
    case 2
    time_arr={'12-12-17';'12-55-11';'12-55-11';'13-01-18';...
              '13-01-18';'13-05-21';'13-05-21';'13-10-11';'13-10-11';...
              '13-14-28';'13-17-00';'13-17-00';'08-24-00';'08-29-36';...
              '08-58-51';'09-51-30';'10-07-23';'10-58-11';'14-40-01';...
              '14-49-10'};
    start_arr=[3,3,64,3,54,3,48,15,80,3,3,64,3,10,3,3,3,3,5,23];
    end_arr=[50,61,127,51,88,45,110,77,110,57,61,118,...
                                               36,35,35,38,43,29,63,66];
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
elseif flight==40030
    
    switch band
    case 1
    time_arr={'13-53-13';'13-53-13';'13-54-36';'14-23-09';'14-23-09';...
              '14-41-29';'14-41-29';'14-45-38';'14-45-38';'14-49-42';...
              '14-49-42';'15-01-10';'15-01-10';'15-10-18';'15-10-18';...
              '15-26-56';'15-29-12';'15-29-12'};
    start_arr=[3,21,3,3,55,3,49,19,68,3,28,3,33,3,51,13,3,66];
    end_arr=[18,34,30,52,79,46,78,65,130,22,83,30,71,48,102,71,63,105];
    case 2
    time_arr={'13-53-13';'13-53-13';'13-54-36';'14-23-09';'14-23-09';...
              '14-41-29';'14-41-29';'14-45-38';'14-45-38';'14-49-42';...
              '15-01-10';'15-01-10';'15-10-18';'15-10-18';'15-26-56';...
              '15-29-12';'15-29-12'};
    start_arr=[3,21,3,3,58,3,50,19,68,3,3,59,13,78,3,3,28];
    end_arr=[18,33,30,55,79,47,78,65,131,23,56,70,75,102,36,25,90];
    end

else
    disp(sprintf('Wrong flight number!'));
end

return 