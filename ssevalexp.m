function sse = ssevalexp(x,xdata,ydata)
a1 = x(1);
% c0 = x(2);
c1 = x(2);
% a2 = x(3);
% c2 = x(4);
% ex = @(x,a,c0,c) c0*ones(size(x))+c.*(1-exp(-x.^2./a^2)).*(x>0)+(0.*(x==0));
% sph = @(x,a,c0,c) (c0*ones(size(x))+c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3))).*(x>0 & x<=a)+(c0+c).*(x>a)+c0.*(x==0);
sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


% sse = nanmean((ydata - sph(xdata,a,c0,c)).^2);
sse = nanmean((ydata -ex(xdata,a1,c1)).^2);
end