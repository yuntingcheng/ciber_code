function dt=get_dark_times(flight,band,field)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gives dark current data time frames for each field(integration time).
%Return full and 1st, 2nd time.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===================================================%
%%%%%%%%%%%%%%%%%%%%%%36265%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================%

if flight==36265
    switch field
    %%%%%%%%%%%%%%%%% 36265 NEP %%%%%%%%%%%%%%%%%%%%    
    case 1
    if band==1
    nfr=35;nfrhalf=17;
    frinfo.name='NEP';frinfo.nfr=35;frinfo.nfrhalf=17;
    frinfo.time{1}='12-35-09';frinfo.start(1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=3;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=20;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=85;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=3;
    elseif band==2
    nfr=35;nfrhalf=17;
    frinfo.name='NEP';frinfo.nfr=35;frinfo.nfrhalf=17;
    frinfo.time{1}='12-35-09';frinfo.start(1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=3;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=19;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=84;
    frinfo.time{end+1}='20-00-01';frinfo.start(end+1)=3;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=3;
    end

    %%%%%%%%%%%%%%%%% 36265 BootesA %%%%%%%%%%%%%%%%%%%%    
    case 2
    if band==1
    nfr=30;nfrhalf=15;
    frinfo.name='BootesA';frinfo.nfr=30;frinfo.nfrhalf=15;
    frinfo.time{1}='12-35-09';frinfo.start(1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=33;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=20;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=50;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=85;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=3;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=33;
    elseif band==2
    nfr=30;nfrhalf=15;
    frinfo.name='BootesA';frinfo.nfr=30;frinfo.nfrhalf=15;
    frinfo.time{1}='12-35-09';frinfo.start(1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=3;
    frinfo.time{end+1}='18-10-47';frinfo.start(end+1)=33;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=19;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=49;
    frinfo.time{end+1}='18-18-01';frinfo.start(end+1)=84;
    frinfo.time{end+1}='20-00-01';frinfo.start(end+1)=3;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=3;
    frinfo.time{end+1}='21-11-35';frinfo.start(end+1)=33;
    end        
    end
end
%===================================================%
%%%%%%%%%%%%%%%%%%%%%%36277%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================%
if flight==36277
    switch field
    %%%%%%%%%%%%%%%%% 36277 SWIRE %%%%%%%%%%%%%%%%%%%%    
    case 1
    if band==1
    nfr=23;nfrhalf=11;
    frinfo.name='SWIRE';frinfo.nfr=23;frinfo.nfrhalf=11;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-12-17';frinfo.start(end+1)=26;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=26;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=65;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=88;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=26;    
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=71;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=16;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=39;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=81;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=26;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=26;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=87;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=8;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=26;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=6;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=29;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23;    
    elseif band==2
    nfr=23;nfrhalf=11;
    frinfo.name='SWIRE';frinfo.nfr=23;frinfo.nfrhalf=11;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-12-17';frinfo.start(end+1)=26;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=26;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=64;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=87;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=26;    
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=71;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=15;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=38;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=80;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=26;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=26;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=87;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=10;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=5;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23; 
    end
    %%%%%%%%%%%%%%%%% 36277 NEP %%%%%%%%%%%%%%%%%%%%    
    case 2
    if band==1
    nfr=24;nfrhalf=12;
    frinfo.name='NEP';frinfo.nfr=24;frinfo.nfrhalf=12;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-12-17';frinfo.start(end+1)=27;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=27;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=65;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=89;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=27;    
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=72;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=16;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=40;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=81;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=27;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=27;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=88;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=8;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=27;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=6;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=30;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23;    
    elseif band==2
    nfr=24;nfrhalf=12;
    frinfo.name='NEP';frinfo.nfr=24;frinfo.nfrhalf=12;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-12-17';frinfo.start(end+1)=27;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=27;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=64;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=88;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=27;    
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=72;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=15;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=39;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=80;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=27;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=27;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=88;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=10;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=5;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=29;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23;    
    end      
    %%%%%%%%%%%%%%%%% 36277 BootesB %%%%%%%%%%%%%%%%%%%%    
    case 3
    if band==1
    nfr=25;nfrhalf=12;
    frinfo.name='BootesB';frinfo.nfr=25;frinfo.nfrhalf=12;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=28;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=65;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=90;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=73;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=16;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=41;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=81;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=28;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=28;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=89;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=8;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-32-16';frinfo.start(end+1)=28;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=6;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=31;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23;    
    elseif band==2
    nfr=25;nfrhalf=12;
    frinfo.name='BootesB';frinfo.nfr=25;frinfo.nfrhalf=12;
    frinfo.time{1}='12-12-17';frinfo.start(1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=28;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=64;
    frinfo.time{end+1}='12-55-11';frinfo.start(end+1)=89;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-01-18';frinfo.start(end+1)=54;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=48;
    frinfo.time{end+1}='13-05-21';frinfo.start(end+1)=73;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=15;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=40;
    frinfo.time{end+1}='13-10-11';frinfo.start(end+1)=80;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-14-28';frinfo.start(end+1)=28;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=28;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=64;
    frinfo.time{end+1}='13-17-00';frinfo.start(end+1)=89;
    frinfo.time{end+1}='08-24-00';frinfo.start(end+1)=3;
    frinfo.time{end+1}='08-29-36';frinfo.start(end+1)=10;
    frinfo.time{end+1}='08-58-51';frinfo.start(end+1)=3;
    frinfo.time{end+1}='09-51-30';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-07-23';frinfo.start(end+1)=3;
    frinfo.time{end+1}='10-58-11';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=5;
    frinfo.time{end+1}='14-40-01';frinfo.start(end+1)=30;
    frinfo.time{end+1}='14-49-10';frinfo.start(end+1)=23;    
    end   
    end
end
%===================================================%
%%%%%%%%%%%%%%%%%%%%%%40030%%%%%%%%%%%%%%%%%%%%%%%%%%
%===================================================%
if flight==40030
    switch field        
    %%%%%%%%%%%%%%%%% 40030 DGL %%%%%%%%%%%%%%%%%%%%
    case 1     
    if band==1
    nfr=38;nfrhalf=19;
    frinfo.name='DGL';frinfo.nfr=38;frinfo.nfrhalf=19;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    
    elseif band==2
    nfr=39;nfrhalf=19;
    frinfo.name='DGL';frinfo.nfr=39;frinfo.nfrhalf=19;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;    
    end
    %%%%%%%%%%%%%%%%% 40030 NEP %%%%%%%%%%%%%%%%%%%%
    case 2     
    if band==1
    nfr=35;nfrhalf=17;
    frinfo.name='NEP';frinfo.nfr=35;frinfo.nfrhalf=17;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=35;nfrhalf=17;
    frinfo.name='NEP';frinfo.nfr=35;frinfo.nfrhalf=17;
    %frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    %frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    %this one goes wrong in the filtering
    frinfo.time{1}='14-45-38';frinfo.start(1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;    
    end
    %%%%%%%%%%%%%%%%% 40030 Lockman %%%%%%%%%%%%%%%%%%%%
    case 3     
    if band==1
    nfr=28;nfrhalf=14;
    frinfo.name='Lockman';frinfo.nfr=28;frinfo.nfrhalf=14;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=96;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=41;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=31;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=29;nfrhalf=14;
    frinfo.name='Lockman';frinfo.nfr=29;frinfo.nfrhalf=14;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=97;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=42;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=57;
    end
    %%%%%%%%%%%%%%%%% 40030 elat10 %%%%%%%%%%%%%%%%%%%%
    case 4     
    if band==1
    nfr=25;nfrhalf=12;
    frinfo.name='elat10';frinfo.nfr=25;frinfo.nfrhalf=12;
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=55;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=49;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=93;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=53;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=76;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=38;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=25;nfrhalf=12;
    frinfo.name='elat10';frinfo.nfr=25;frinfo.nfrhalf=12;        
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=93;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=38;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=78;    
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=53;
    end
    %%%%%%%%%%%%%%%%% 40030 elat30 %%%%%%%%%%%%%%%%%%%%
    case 5    
    if band==1 
    nfr=26;nfrhalf=13;
    frinfo.name='elat30';frinfo.nfr=26;frinfo.nfrhalf=13;
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=49;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=94;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=54;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=77;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=39;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=29;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=26;nfrhalf=13;
    frinfo.name='elat30';frinfo.nfr=26;frinfo.nfrhalf=13;
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=29;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=94;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=29;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=39;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=54;
    end
    %%%%%%%%%%%%%%%%% 40030 BootesB %%%%%%%%%%%%%%%%%%%%
    case 6 
    if band==1
    nfr=30;nfrhalf=15;
    frinfo.name='BootesB';frinfo.nfr=30;frinfo.nfrhalf=15;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=49;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=98;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=41;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=29;nfrhalf=14;
    frinfo.name='BootesB';frinfo.nfr=29;frinfo.nfrhalf=14;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=97;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=42;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=57;    
    end
    %%%%%%%%%%%%%%%%% 40030 BootesA %%%%%%%%%%%%%%%%%%%%
    case 7     
    if band==1
    nfr=29;nfrhalf=14;
    frinfo.name='BootesA';frinfo.nfr=29;frinfo.nfrhalf=14;
    frinfo.time{1}='14-23-09';frinfo.start(1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=49;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=97;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=42;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=32;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;        
    elseif band==2
    nfr=28;nfrhalf=14;
    frinfo.name='BootesA';frinfo.nfr=28;frinfo.nfrhalf=14;
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=96;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=41;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=56;
    end    
    %%%%%%%%%%%%%%%%% 40030 SWIRE %%%%%%%%%%%%%%%%%%%%
    case 8 
    if band==1
    nfr=26;nfrhalf=13;
    frinfo.name='SWIRE';frinfo.nfr=26;frinfo.nfrhalf=13;
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=49;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=94;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=28;
    frinfo.time{end+1}='14-49-42';frinfo.start(end+1)=54;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=33;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=51;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=77;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=39;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=29;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=66;
    elseif band==2
    nfr=26;nfrhalf=13;
    frinfo.name='SWIRE';frinfo.nfr=26;frinfo.nfrhalf=13;        
    frinfo.time{1}='13-54-36';frinfo.start(1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-23-09';frinfo.start(end+1)=29;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=3;
    frinfo.time{end+1}='14-41-29';frinfo.start(end+1)=50;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=19;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=68;
    frinfo.time{end+1}='14-45-38';frinfo.start(end+1)=94;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-01-10';frinfo.start(end+1)=29;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=13;
    frinfo.time{end+1}='15-10-18';frinfo.start(end+1)=39;
    frinfo.time{end+1}='15-26-56';frinfo.start(end+1)=3;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=28;
    frinfo.time{end+1}='15-29-12';frinfo.start(end+1)=54;
    end    
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt.name=frinfo.name;dt.nfr=frinfo.nfr;dt.nfrhalf=frinfo.nfrhalf;
    dt.time{1}=frinfo.time{1};
    dt.frdown=frinfo.start(1);
    dt.frup=dt.frdown+nfr-1;
    dt.frdown1=dt.frdown;
    dt.frup1=dt.frdown1+nfrhalf-1;
    dt.frdown2=dt.frup1+1;
    dt.frup2=dt.frdown2+nfrhalf-1;
    
    for i=2:numel(frinfo.time)
    dt.time{end+1}=frinfo.time{i};
    dt.frdown(end+1)=frinfo.start(i);
    dt.frup(end+1)=dt.frdown(end)+nfr-1;
    dt.frdown1(end+1)=dt.frdown(end);
    dt.frup1(end+1)=dt.frdown1(end)+nfrhalf-1;
    dt.frdown2(end+1)=dt.frup1(end)+1;
    dt.frup2(end+1)=dt.frdown2(end)+nfrhalf-1;
    end       
return
