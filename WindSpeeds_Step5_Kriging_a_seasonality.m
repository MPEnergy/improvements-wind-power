%% Kriging of the alpha1, alpha2 AR(2) parameters

% Read the parameters
clear all
p = xlsread('NewSeasonality.xlsx', 'Yearly');


lats = p(:,2);
lons = p(:,3);
A0 = p(:,4);
A1 = p(:,5);
A2 = p(:,6);
A3 = p(:,7);
A4 = p(:,8);
A5 = p(:,9);
A6 = p(:,10);
A7 = p(:,11);
A8 = p(:,12);
A9 = p(:,13);
A10 = p(:,14);
A11 = p(:,15);
A12 = p(:,16);
%% Make the semivariogram
for i= 1:85
    for j = 1:85
        [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
    end 
end

Ns = unique(D);
edges = 0:11:800;
[Nb, BIN] = histc(Ns, edges);


%% A0
[g0, cf0] = calculate_semivariogram(A0, edges, D);

%% A1
[g1, cf1] = calculate_semivariogram(A1, edges, D);

%% A2
[g2, cf2] = calculate_semivariogram(A2, edges, D);

%% A3
[g3, cf3] = calculate_semivariogram(A3, edges, D);

%% A4
[g4, cf4] = calculate_semivariogram(A4, edges, D);

%% A5
[g5, cf5] = calculate_semivariogram(A5, edges, D);

%% A6
[g6, cf6] = calculate_semivariogram(A6, edges, D);
%% A7
[g7, cf7] = calculate_semivariogram(A7, edges, D);

%% A8
[g8, cf8] = calculate_semivariogram(A8, edges, D);

%% A9
[g9, cf9] = calculate_semivariogram(A9, edges, D);

%% A10
[g10, cf10] = calculate_semivariogram(A10, edges, D);

%% A11
[g11, cf11] = calculate_semivariogram(A11, edges, D);

%% A12
[g12, cf12] = calculate_semivariogram(A12, edges, D);


%% Plot
% for i=0:12
%     figure(i+1);
%     field = sprintf('g%d',i);
%     scatter(edges, eval(field))
%     xlabel('Distance(km)')
%     ylabel('Semivariogram')
%     title(sprintf('Semivariogram a_{%d}', i))
% end


%% Fit the spatial averaged semivariogram
% load('Semivariograms_a_coeffs.mat','edges','g0','g1','g2','g3','g4','g5','g6','g7','g8','g9','g10','g11', 'g12','cf0','cf1','cf2','cf3','cf4','cf5','cf6','cf7','cf8','cf9','cf10','cf11', 'cf12')
% A5
xdata = edges';
ydata = g5';
fun1 = @(x)ssevalexp1(x,xdata,ydata);


%x0 = [400, 100]; %% [400,10, 50, -20]
x0 = [600,0.0004, 200, -5];
options = optimset('MaxFunEvals',100000000);
bestx1 = fminsearch(fun1,x0, options);
% sph = @(x,a,c0,c) c0*ones(size(x))+c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+(c0+c).*(x>a)+c0.*(x==0);
sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);

%Check the fit
a1 = bestx1(1);
c1 = bestx1(2);
a2 = bestx1(3);
c2 = bestx1(4);


yfit = sph(xdata, a1, c1)+gaus(xdata,a2,c2);
figure;
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram (m/s)')
title(sprintf('Sine Hole Effect and Exponential function'))
legend('a_0 semivariogram','Fitted semivariogram')
hold off

sqrt(nanmean((ydata-yfit).^2))

%% a0
a1 = 285.0907;
c1 = 0.3950;
a2 = 196.0287;
c2 = -0.2474;

xdata = edges';
ydata = g0';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fita0 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+sinehole(h,a2,c2));  %the fitted semivariogram function
yfit = fita0(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_0'))
legend('a_0 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf0';
fitcov0 = @(h,a1,c1,a2,c2) (0.1198 -fita0(h,a1,c1,a2,c2));

yfit = 0.1198 -fita0(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_0'))
legend('a_0 covariance function','Fitted covariance')
hold off
nugget0 = 0.1198;

%% a1
a1 = 927.6168;
c1 = -1.1021;
a2 = 1190.3;
c2 = 2.3935;

xdata = edges';
ydata = g1';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita1 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fita1(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_1'))
legend('a_1 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf1';
fitcov1 = @(h,a1,c1,a2,c2) (0.01615 -fita1(h,a1,c1,a2,c2));

yfit = 0.01615 -fita1(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_1'))
legend('a_1 covariance function','Fitted covariance')
hold off
nugget1 = 0.01615; 
%% a2
a1 = 715.5994;
c1 = -0.00065932;
a2 = 55.1084;
c2 = 0.0014;

xdata = edges';
ydata = g2';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita2 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fita2(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_2'))
legend('a_2 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf2';
fitcov2 = @(h,a1,c1,a2,c2) (0.0009625-fita2(h,a1,c1,a2,c2));

yfit = 0.0009625-fita2(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_2'))
legend('a_2 covariance function','Fitted covariance')
hold off
nugget2 = 0.0009625;

%% a3
a1 = 98.3195;
c1 = 0.0010;
a2 = 814.5254;
c2 = 0.00040862;

xdata = edges';
ydata = g3';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita3 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fita3(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_3'))
legend('a_3 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf3';
fitcov3 = @(h,a1,c1,a2,c2) (0.00109-fita3(h,a1,c1,a2,c2));

yfit = 0.00109-fita3(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_3'))
legend('a_3 covariance function','Fitted covariance')
hold off
nugget3 = 0.00109;

%% a4
a1 = 564.5545;
c1 = 0.0027;
a2 = 407.8630;
c2 = -0.0018;

xdata = edges';
ydata = g4';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita4 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+sinehole(h,a2,c2));  %the fitted semivariogram function
yfit = fita4(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_4'))
legend('a_4 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf4';
fitcov4 = @(h,a1,c1,a2,c2) (0.000623-fita4(h,a1,c1,a2,c2));

yfit = 0.000623-fita4(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_4'))
legend('a_4 covariance function','Fitted covariance')
hold off
nugget4 = 0.000623;

%% a5
a1 = 511.6893;
c1 = 0.0014;
a2 = 231.3095;
c2 = -0.0011;

xdata = edges';
ydata = g5';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita5 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+gaus(h,a2,c2));  %the fitted semivariogram function
yfit = fita5(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_5'))
legend('a_5 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf5';
fitcov5 = @(h,a1,c1,a2,c2) (0.000254-fita5(h,a1,c1,a2,c2));

yfit = 0.000254-fita5(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_5'))
legend('a_5 covariance function','Fitted covariance')
hold off
nugget5 = 0.000254;

%% a6
a1 = 156.6367;
c1 = 0.00063892;
a2 = 291.3089;
c2 = -0.00032461;


xdata = edges';
ydata = g6';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita6 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sinehole(h,a2,c2));  %the fitted semivariogram function
yfit = fita6(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_6'))
legend('a_6 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf6';
fitcov6 = @(h,a1,c1,a2,c2) (0.0002613-fita6(h,a1,c1,a2,c2));

yfit = 0.0002613-fita6(xdata,a1,c1,a2,c2);
subplot(2,1,2);
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_6'))
legend('a_6 covariance function','Fitted covariance')
hold off
nugget6 = 0.0002613;

%% a7
a1 = 98.3912;
c1 = 0.00019393;
a2 = 4828.4;
c2 = 0.01;


xdata = edges';
ydata = g7';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita7 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+gaus(h,a2,c2));  %the fitted semivariogram function
yfit = fita7(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_7'))
legend('a_7 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf7';
fitcov7 = @(h,a1,c1,a2,c2) (0.00022-fita7(h,a1,c1,a2,c2));

yfit = 0.00022-fita7(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_7'))
legend('a_7 covariance function','Fitted covariance')
hold off
nugget7 = 0.00022;
%% a8
a1 = 208.3804;
c1 = 0.00029236;
a2 = 249.7290;
c2 = -0.0001579;


xdata = edges';
ydata = g8';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita8 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+gaus(h,a2,c2));  %the fitted semivariogram function
yfit = fita8(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_8'))
legend('a_8 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf8';
fitcov8 = @(h,a1,c1,a2,c2) (0.000173-fita8(h,a1,c1,a2,c2));
yfit = 0.000173-fita8(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_8'))
legend('a_8 covariance function','Fitted covariance')
hold off
nugget8 = 0.000173;

%% a9
a1 = 570.0844;
c1 = 0.00030183;
a2 = 357190;
c2 = -0.0468;


xdata = edges';
ydata = g9';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fita9 = @(h,a1,c1,a2,c2) (pentasph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fita9(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_9'))
legend('a_9 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf9';
fitcov9 = @(h,a1,c1,a2,c2) (0.0001904-fita9(h,a1,c1,a2,c2)); 
yfit = 0.0001904-fita9(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_9'))
legend('a_9 covariance function','Fitted covariance')
hold off
nugget9 = 0.0001904;

%% a10
a1 = 554.0978;
c1 = 0.00051255;
a2 = 4354.8;
c2 = -0.0032;


xdata = edges';
ydata = g10';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita10 = @(h,a1,c1,a2,c2) (pentasph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fita10(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_{10}'))
legend('a_{10} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf10';
fitcov10 = @(h,a1,c1,a2,c2) (0.000166-fita10(h,a1,c1,a2,c2)); 
yfit = 0.000166-fita10(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_{10}'))
legend('a_{10} covariance function','Fitted covariance')
hold off
nugget10 = 0.000166;

%% a11
a1 = 556.6596;
c1 = 0.0011;
a2 = 419.6655;
c2 = -0.00072825;


xdata = edges';
ydata = g11';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita11 = @(h,a1,c1,a2,c2)(sph(h,a1,c1)+sinehole(h,a2,c2));  %the fitted semivariogram function
yfit = fita11(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_{11}'))
legend('a_{11} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf11';
fitcov11 = @(h,a1,c1,a2,c2) (0.0002626-fita11(h,a1,c1,a2,c2)); 
yfit = 0.0002626-fita11(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_{11}'))
legend('a_{11} covariance function','Fitted covariance')
hold off
nugget11 = 0.0002626;

%% a12
a1 = 629.1118;
c1 = 0.00083232;
a2 = 322.9730;
c2 = -0.00056099;


xdata = edges';
ydata = g12';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0


fita12 = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+gaus(h,a2,c2));  %the fitted semivariogram function
yfit = fita12(xdata,a1,c1,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title(sprintf('Semivariogram of a_{12}'))
legend('a_{12} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cf12';

fitcov12 = @(h,a1,c1,a2,c2) (0.0002073-fita12(h,a1,c1,a2,c2)); 
yfit = 0.0002073-fita12(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
title(sprintf('Covariance function of a_{12}'))
legend('a_{12} covariance function','Fitted covariance')
hold off
nugget12 = 0.0002073;


%% Kriging for a0
delta = 0.01;
newlats = min(lats):0.01:max(lats);
newlons = min(lons):0.01:max(lons);
[lat, lon] = meshgrid(newlats, newlons);


newlat = [];
newlon = [];
for r=1:751
    for c=1:301
        if lat(r,c)~=0 && lon(r,c)~=0
            newlat = [newlat lat(r,c)];
            newlon = [newlon lon(r,c)];
        end
    end
end
newlat = newlat';
newlon = newlon';


%Limit to the right polygon


idx1 = newlon<-newlat-84.5 & newlat<33.5 & newlat>=32.5;


idx2 = newlon<-2*newlat-51 & newlat>=33.5 & newlat<34;


idx3 = newlon<-3*newlat-17 & newlat>=34 & newlat<34.5;


idx4 = newlon<-120.5 & newlat>=34.5  & newlat<35;


idx5 = newlon<-newlat-85.5 & newlat>=35 & newlat<=35.5;

newlon = newlon(~(idx1 | idx2 | idx3 | idx4 | idx5));
newlat = newlat(~(idx1 | idx2 | idx3 | idx4 | idx5));





% Kriging per se
%% a0
n = length(newlat);
krigged_a0 = zeros(n,1);

a1 = 285.0907;
c1 = 0.3950;
a2 = 196.0287;
c2 = -0.2474;

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A0, fitcov0, fita0, nugget0, sph, sinehole, a1, c1, a2, c2);
end


krigged_a0 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a0, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_0')

save('KriggedResults_a0.mat','A0','krigged_a0', 'lat','lon')

%% a1
krigged_a1 = zeros(n,1);

a1 = 927.6168;
c1 = -1.1021;
a2 = 1190.3;
c2 = 2.3935;

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A1, fitcov1, fita1, nugget1, sph, ex, a1, c1, a2, c2);
end
krigged_a1 = krigged_a';


figure
scatter(newlon,newlat, 20, krigged_a1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_1')

save('KriggedResults_a1.mat','A1','krigged_a1', 'lat','lon')

%% a2
krigged_a2 = zeros(n,1);

a1 = 715.5994;
c1 = -0.00065932;
a2 = 55.1084;
c2 = 0.0014;

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A2, fitcov2, fita2, nugget2, sph, ex, a1, c1, a2, c2);
end
krigged_a2 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a2, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_2')

save('KriggedResults_a2.mat','A2','krigged_a2', 'lat','lon')

n = length(newlat);
%% a3
a1 = 98.3195;
c1 = 0.0010;
a2 = 814.5254;
c2 = 0.00040862;

krigged_a3 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A3, fitcov3, fita3, nugget3, sph, ex, a1, c1, a2, c2);
end
krigged_a3 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a3, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_3')

save('KriggedResults_a3.mat','A3','krigged_a3', 'lat','lon')

%% a4
a1 = 564.5545;
c1 = 0.0027;
a2 = 407.8630;
c2 = -0.0018;


krigged_a4 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A4, fitcov4, fita4, nugget4, sph, sinehole, a1, c1, a2, c2);
end
krigged_a4 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a4, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_4')

save('KriggedResults_a4.mat','A4','krigged_a4', 'lat','lon')


%% a5
a1 = 511.6893;
c1 = 0.0014;
a2 = 231.3095;
c2 = -0.0011;

krigged_a5 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A5, fitcov5, fita5, nugget5, sph, gaus, a1, c1, a2, c2);
end
krigged_a5 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a5, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_5')

save('KriggedResults_a5.mat','A5','krigged_a5', 'lat','lon')


%% a6
a1 = 156.6367;
c1 = 0.00063892;
a2 = 291.3089;
c2 = -0.00032461;

krigged_a6 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A6, fitcov6, fita6, nugget6, gaus, sinehole, a1, c1, a2, c2);
end
krigged_a6 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a6, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_6')

save('KriggedResults_a6.mat','A6','krigged_a6', 'lat','lon')


%% a7
a1 = 98.3912;
c1 = 0.00019393;
a2 = 4828.4;
c2 = 0.01;

krigged_a7 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A7, fitcov7, fita7, nugget7, sph, gaus, a1, c1, a2, c2);
end
krigged_a7 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a7, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_7')

save('KriggedResults_a7.mat','A7','krigged_a7', 'lat','lon')

%% a8
a1 = 208.3804;
c1 = 0.00029236;
a2 = 249.7290;
c2 = -0.0001579;

krigged_a8 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A8, fitcov8, fita8, nugget8, sph, gaus, a1, c1, a2, c2);
end
krigged_a8 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a8, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_8')

save('KriggedResults_a8.mat','A8','krigged_a8', 'lat','lon')


%% a9
a1 = 570.0844;
c1 = 0.00030183;
a2 = 357190;
c2 = -0.0468;

krigged_a9 = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A9, fitcov9, fita9, nugget9, pentasph, ex, a1, c1, a2, c2);
end
krigged_a9 = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_a9, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_9')

save('KriggedResults_a9.mat','A9','krigged_a9', 'lat','lon')

%% a10
a1 = 554.0978;
c1 = 0.00051255;
a2 = 4354.8;
c2 = -0.0032;

krigged_adiez = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A10, fitcov10, fita10, nugget10, pentasph, ex, a1, c1, a2, c2);
end
krigged_adiez = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_adiez, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_{10}')

save('KriggedResults_a10.mat','A10','krigged_adiez', 'lat','lon')


%% a11
a1 = 556.6596;
c1 = 0.0011;
a2 = 419.6655;
c2 = -0.00072825;


krigged_aonce = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A11, fitcov11, fita11, nugget11, sph, sinehole, a1, c1, a2, c2);
end
krigged_aonce = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_aonce, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_{11}')

save('KriggedResults_a11.mat','A11','krigged_aonce', 'lat','lon')

%% a12
a1 = 629.1118;
c1 = 0.00083232;
a2 = 322.9730;
c2 = -0.00056099;


krigged_adoce = zeros(n,1);

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, A12, fitcov12, fita12, nugget12, sph, gaus, a1, c1, a2, c2);
end
krigged_adoce = krigged_a';

figure
scatter(newlon,newlat, 20, krigged_adoce, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_{12}')

save('KriggedResults_a12.mat','A12','krigged_adoce', 'lat','lon')

%%
figure
scatter(newlon,newlat, 20, krigged_a2, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Krigging results for a_2')

figure
scatter(lons,lats, 60, a0, 'filled');
colormap(hsv) %can be removed
xlabel('Lat')
ylabel('Lon')
colorbar
title('a_0')

%%

for p = 1:n
    krigged_a(p) = kriging(p,newlat, newlon, lats, lons, a0, fitcov0, fita0, gaus, sinehole, a1, c1, a2, c2);
end


kirgged_a0 = krigged_a;


function krigged_alpha = local_kriging(p, newlat, newlon, lats, lons, alpha, covfun, fitalpha, fun1, fun2, param1, param2, param3, param4)
if newlat(p)-round(newlat(p))~=0.5 || newlon(p)-round(newlon(p))~=0.5
    %C*lambda=v
    n = length(newlat);
    v = [];
    for i=1:n
        [dist,~] = lldistkm([newlat(i) newlon(i)], [newlat(p) newlon(p)]); 
       v = [v covfun(dist, param1, param2, param3, param4)]; 
    end
    
    
    
    
    
    
else %it is already a point in the grid
    idxlat = find(lats==newlat(p));
    idxlon = find(lons==newlon(p));
    if idxlat==idxlon
        krigged_alpha = alpha1(idxlat);
    else
        sprintf('Wrong lat %d, lon %d', newlat(p), newlon(p));
        return;
    end
end
end

function [g, cf] = calculate_semivariogram(z, edges, D)
    for k=1:length(edges)-1
        [idxrow, idxcol] = find(D>=edges(k) & D<edges(k+1));

        gamma = 0;
        for i=1:length(idxrow)
            gamma = gamma + (z(idxrow(i))-z(idxcol(i)))^2;
        end
        num = length(idxrow);
        g(k) = gamma/(2*num);
        C = cov(z(idxrow),z(idxcol));
        cf(k) = C(1,2);
    end
    k = length(edges);
    [idxrow, idxcol] = find(D>=edges(k-1));
    gamma = 0;
    for i=1:length(idxrow)
        gamma = gamma + (z(idxrow(i))-z(idxcol(i)))^2
    end
    C = cov(z(idxrow),z(idxcol));
    cf(k) = C(1,2);
    num = length(idxrow);
    g(k) = gamma/(2*num);
end



