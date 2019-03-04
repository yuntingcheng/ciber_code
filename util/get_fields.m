function fields=get_fields(flight,inst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Gives flight number(36265,36277,40030) and band (1 or 2), 
%gives field name, Nfr, and Nfr_half in flight data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch flight
case 36265
    fields(1).name='NEP';fields(1).nfr=35;fields(1).nfrhalf=17;
    fields(2).name='BootesA';fields(2).nfr=30;fields(2).nfrhalf=15;
case 36277
    fields(1).name='SWIRE';fields(1).nfr=23;fields(1).nfrhalf=11;
    fields(2).name='NEP';fields(2).nfr=24;fields(2).nfrhalf=12;
    fields(3).name='BootesB';fields(3).nfr=25;fields(3).nfrhalf=12;
case 40030
    if inst==1
    fields(1).name='DGL';fields(1).nfr=38;fields(1).nfrhalf=19;
    fields(2).name='NEP';fields(2).nfr=35;fields(2).nfrhalf=17;
    fields(3).name='Lockman';fields(3).nfr=28;fields(3).nfrhalf=14;
    fields(4).name='elat10';fields(4).nfr=25;fields(4).nfrhalf=12;
    fields(5).name='elat30';fields(5).nfr=26;fields(5).nfrhalf=13;
    fields(6).name='BootesB';fields(6).nfr=30;fields(6).nfrhalf=15;
    fields(7).name='BootesA';fields(7).nfr=29;fields(7).nfrhalf=14;
    fields(8).name='SWIRE';fields(8).nfr=26;fields(8).nfrhalf=13;
    elseif inst==2
    fields(1).name='DGL';fields(1).nfr=39;fields(1).nfrhalf=19;
    fields(2).name='NEP';fields(2).nfr=35;fields(2).nfrhalf=17;
    fields(3).name='Lockman';fields(3).nfr=29;fields(3).nfrhalf=14;
    fields(4).name='elat10';fields(4).nfr=25;fields(4).nfrhalf=12;
    fields(5).name='elat30';fields(5).nfr=26;fields(5).nfrhalf=13;
    fields(6).name='BootesB';fields(6).nfr=29;fields(6).nfrhalf=14;
    fields(7).name='BootesA';fields(7).nfr=28;fields(7).nfrhalf=14;
    fields(8).name='SWIRE';fields(8).nfr=26;fields(8).nfrhalf=13;
    end  
end
end
