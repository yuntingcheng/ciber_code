function sse=my2dgaussfit(params,Input,Actual_Output)

sig = params(1);
x0 = params(2);
y0 = params(3);
A = params(4);

s = size(Actual_Output);

x=zeros(s(1));
y=x;

for i=1:s(1)
    x(i,:)=i;
    y(:,i)=i;
end

x=x(:);
y=y(:);
Actual_Output = Actual_Output(:);
Fitted_Curve = A*exp(- ( ((x - x0).^2)/(2*sig^2) + ((y -y0).^2)/(2*sig^2)  )  );
Fitted_Curve = Fitted_Curve(:);

Error_Vector=Fitted_Curve - Actual_Output;
% When curvefitting, a typical quantity to
% minimize is the sum of squares error
sse=sum(Error_Vector.^2);
% You could also write sse as
% sse=Error_Vector(:)'*Error_Vector(:);

return