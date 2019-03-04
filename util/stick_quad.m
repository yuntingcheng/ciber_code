function map = stick_quad(quadAname)

map = zeros(1024);
quad_arr=['A','B','C','D'];
for iquad=1:4
    quad=quad_arr(iquad);
    quadname = strrep(quadAname,'_A_',strcat('_',quad,'_'));
    qmap=fits_read(quadname);

    if iquad==1
        map(1:512,1:512)=qmap;
    elseif iquad==2
        map(513:1024,1:512)=qmap;
    elseif iquad==3
        map(1:512,513:1024)=qmap;
    else
        map(513:1024,513:1024)=qmap;
    end
end


return