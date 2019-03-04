function sse=gaussfit2d(params,data)

sigx = params(1);
sigy = params(2);
x0 = params(3);
y0 = params(4);
A = params(5);

s = size(data);

x=zeros(s(1));
y=x;

for i=1:s(1)
    x(i,:)=i;
    y(:,i)=i;
end

x=x(:);
y=y(:);
Fitted_Curve = A*exp(-(((x - x0).^2)./(2*sigx^2) + ((y -y0).^2)./(2*sigy^2)));
Error_Vector=Fitted_Curve(:) - data(:);
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);

return
