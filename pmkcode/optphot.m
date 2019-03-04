function [flux,neff] = optphot(source,beam)

beam=beam./sum(beam(:));


weight = beam./sum((beam(:).^2)) ;
flux = sum(weight(:).*source(:));


s = sum(beam(:));
x = beam./s;

w2=sum(x(:).^2);

neff = 1./w2;

end