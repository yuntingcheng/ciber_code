function noisey = awgn_pk(x,sig)

noisey = x.*0;

n = awgn(noisey,1);

st = std(n(:));

n=sig*(n./st);

noisey(:) = x(:)+n(:);

end