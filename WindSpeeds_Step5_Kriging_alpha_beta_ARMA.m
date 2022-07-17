%% Kriging of the alpha1, alpha2, beta1, ...beta24 ARMA(2,24) parameters

% Read the parameters
clear all
p = xlsread('ARMA.xlsx','ARMA(2,24)');

lats = p(:,2);
lons = p(:,3);
alpha1 = p(:,6);
alpha2 = p(:,7);
beta1 = p(:,8);
beta2 = p(:,9);
beta3 = p(:,10);
beta4 = p(:,11);
beta5 = p(:,12);
beta6 = p(:,13);
beta7 = p(:,14);
beta8 = p(:,15);
beta9 = p(:,16);
beta10 = p(:,17);
beta11 = p(:,18);
beta12 = p(:,19);
beta13 = p(:,20);
beta14 = p(:,21);
beta15 = p(:,22);
beta16 = p(:,23);
beta17 = p(:,24);
beta18 = p(:,25);
beta19 = p(:,26);
beta20 = p(:,27);
beta21 = p(:,28);
beta22 = p(:,29);
beta23 = p(:,30);
beta24 = p(:,31);
intercept = p(:,32);

%% Make the semivariogram
for i= 1:85
    for j = 1:85
        [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
    end 
end

Ns = unique(D);
edges = 1:11:800;
[Nb, BIN] = histc(Ns, edges);


%% alpha1
[galpha1, cfalpha1] = calculate_semivariogram(alpha1, edges, D);
%% alpha2
[galpha2, cfalpha2] = calculate_semivariogram(alpha2, edges, D);
%% beta1
[gbeta1, cfbeta1] = calculate_semivariogram(beta1, edges, D);
%% beta2
[gbeta2, cfbeta2] = calculate_semivariogram(beta2, edges, D);
%% beta3
[gbeta3, cfbeta3] = calculate_semivariogram(beta3, edges, D);
%% beta4
[gbeta4, cfbeta4] = calculate_semivariogram(beta4, edges, D);
%% beta5
[gbeta5, cfbeta5] = calculate_semivariogram(beta5, edges, D);
%% beta6
[gbeta6, cfbeta6] = calculate_semivariogram(beta6, edges, D);
%% beta7
[gbeta7, cfbeta7] = calculate_semivariogram(beta7, edges, D);
%% beta8
[gbeta8, cfbeta8] = calculate_semivariogram(beta8, edges, D);
%% beta9
[gbeta9, cfbeta9] = calculate_semivariogram(beta9, edges, D);
%% beta10
[gbeta10, cfbeta10] = calculate_semivariogram(beta10, edges, D);
%% beta11
[gbeta11, cfbeta11] = calculate_semivariogram(beta11, edges, D);
%% beta12
[gbeta12, cfbeta12] = calculate_semivariogram(beta12, edges, D);
%% beta13
[gbeta13, cfbeta13] = calculate_semivariogram(beta13, edges, D);
%% beta14
[gbeta14, cfbeta14] = calculate_semivariogram(beta14, edges, D);
%% beta15
[gbeta15, cfbeta15] = calculate_semivariogram(beta15, edges, D);
%% beta16
[gbeta16, cfbeta16] = calculate_semivariogram(beta16, edges, D);
%% beta17
[gbeta17, cfbeta17] = calculate_semivariogram(beta17, edges, D);
%% beta18
[gbeta18, cfbeta18] = calculate_semivariogram(beta18, edges, D);
%% beta19
[gbeta19, cfbeta19] = calculate_semivariogram(beta19, edges, D);
%% beta20
[gbeta20, cfbeta20] = calculate_semivariogram(beta20, edges, D);
%% beta21
[gbeta21, cfbeta21] = calculate_semivariogram(beta21, edges, D);
%% beta22
[gbeta22, cfbeta22] = calculate_semivariogram(beta22, edges, D);
%% beta23
[gbeta23, cfbeta23] = calculate_semivariogram(beta23, edges, D);
%% beta24
[gbeta24, cfbeta24] = calculate_semivariogram(beta24, edges, D);


% save('Spatial_semivariogram_alphasbetas.mat','galpha1','g2','alpha1','alpha2','lats','lons','edges')

%% Fit the spatial averaged semivariogram

%% Beta24
xdata = edges';
ydata = galpha1';
fun1 = @(x)ssevalexp1(x,xdata,ydata);

% x0 = [100, 0.0004];
 x0 = [1000, -0.1625, 21, 0.35];
 %x0 = [100, 10, 40, 0.089];
options = optimset('MaxFunEvals',10000000000);
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

yfit = pentasph(xdata,a1,c1)+ex(xdata,a2,c2);
figure;
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram (m/s)')
title('Fit of spatial variation for $\alpha_1$', 'Interpreter', 'latex')
legend('Alpha_1 semivariogram','Fitted semivariogram')
hold off

sqrt(nanmean((ydata-yfit).^2))




%% alpha1
a1 = 1003.9;
c1 = -0.1625;
a2 = 21.3262;
c2 = 0.3516;

xdata = edges';
ydata = galpha1';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitalpha1 = @(h,a1,c1,a2,c2) (pentasph(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitalpha1(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\alpha_1$', 'Interpreter', 'latex')
legend('\alpha_1 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfalpha1';
fitcovalpha1 = @(h,a1,c1,a2,c2) (0.2775 -fitalpha1(h,a1,c1,a2,c2));

yfit = 0.2775 -fitalpha1(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\alpha_1 covariance function','Fitted covariance')
hold off
nuggetalpha1 = 0.2775;


%% alpha2
a1 = 350.9542;
c1 = -0.1151;
a2 = 25.6248;
c2 = 0.2793;

xdata = edges';
ydata = galpha2';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitalpha2 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitalpha2(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\alpha_2$', 'Interpreter', 'latex')
legend('\alpha_2 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfalpha2';
fitcovalpha2 = @(h,a1,c1,a2,c2) (0.2282 -fitalpha2(h,a1,c1,a2,c2));

yfit = 0.2282 -fitalpha2(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\alpha_2 covariance function','Fitted covariance')
hold off
nuggetalpha2 = 0.2282;


%% beta1
a1 = 350.4368;
c1 = -0.1371;
a2 = 23.8057;
c2 = 0.3457;

xdata = edges';
ydata = gbeta1';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta1 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta1(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_1$', 'Interpreter', 'latex')
legend('\beta_1 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta1';
fitcovbeta1 = @(h,a1,c1,a2,c2) (0.2853 -fitbeta1(h,a1,c1,a2,c2));

yfit = 0.2853 -fitbeta1(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_1 covariance function','Fitted covariance')
hold off
nuggetbeta1 = 0.2853;


%% beta2
a1 = 256.6674;
c1 = -0.0025;
a2 = 4.9492;
c2 = 0.0121;

xdata = edges';
ydata = gbeta2';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta2 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta2(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_2$', 'Interpreter', 'latex')
legend('\beta_2 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta2';
fitcovbeta2 = @(h,a1,c1,a2,c2) (0.01 -fitbeta2(h,a1,c1,a2,c2));

yfit = 0.01 -fitbeta2(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_2 covariance function','Fitted covariance')
hold off
nuggetbeta2 = 0.01;



%% beta3
a1 = 124.3028;
c1 = -0.0021;
a2 = 32.8914;
c2 = 0.0061;

xdata = edges';
ydata = gbeta3';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta3 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta3(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_3$', 'Interpreter', 'latex')
legend('\beta_3 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta3';
fitcovbeta3 = @(h,a1,c1,a2,c2) (0.004189 -fitbeta3(h,a1,c1,a2,c2));

yfit = 0.004189 -fitbeta3(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_3 covariance function','Fitted covariance')
hold off
nuggetbeta3 = 0.004189;



%% beta4
a1 = 98.1568;
c1 = -0.00043031;
a2 = 19.7806;
c2 = 0.0032;

xdata = edges';
ydata = gbeta4';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta4 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta4(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_4$', 'Interpreter', 'latex')
legend('\beta_4 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta4';
fitcovbeta4 = @(h,a1,c1,a2,c2) (0.00276 -fitbeta4(h,a1,c1,a2,c2));

yfit = 0.00276 -fitbeta4(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_4 covariance function','Fitted covariance')
hold off
nuggetbeta4 = 0.00276;


%% beta5
a1 = 226.4135;
c1 = -0.0051;
a2 = 248.5329;
c2 = 0.0085;

xdata = edges';
ydata = gbeta5';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta5 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta5(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_5$', 'Interpreter', 'latex')
legend('\beta_5 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta5';
fitcovbeta5 = @(h,a1,c1,a2,c2) (0.001989 -fitbeta5(h,a1,c1,a2,c2));

yfit = 0.001989 -fitbeta5(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_5 covariance function','Fitted covariance')
hold off
nuggetbeta5 = 0.001989;


%% beta6
a1 = 180.9366;
c1 = -0.0039;
a2 = 166.5990;
c2 = 0.0058;

xdata = edges';
ydata = gbeta6';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta6 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta6(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_6$', 'Interpreter', 'latex')
legend('\beta_6 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta6';
fitcovbeta6 = @(h,a1,c1,a2,c2) (0.001424 -fitbeta6(h,a1,c1,a2,c2));

yfit = 0.001424 -fitbeta6(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_6 covariance function','Fitted covariance')
hold off
nuggetbeta6 = 0.001424;


%% beta7
a1 = -16.4444;
c1 = 0.0010;
a2 = 397.6020;
c2 = -0.00003242;

xdata = edges';
ydata = gbeta7';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta7 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta7(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_7$', 'Interpreter', 'latex')
legend('\beta_7 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta7';
fitcovbeta7 = @(h,a1,c1,a2,c2) (0.001025 -fitbeta7(h,a1,c1,a2,c2));

yfit = 0.001025 -fitbeta7(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_7 covariance function','Fitted covariance')
hold off
nuggetbeta7 = 0.001025;


%% beta8
a1 = 30.2180;%37.1070;
c1 = 0.00070866;%0.00086581;
a2 = 226.6858;%236.2554;
c2 = -0.000034576;%-0.0001985;

xdata = edges';
ydata = gbeta8';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta8 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta8(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_8$', 'Interpreter', 'latex')
legend('\beta_8 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta8';
fitcovbeta8 = @(h,a1,c1,a2,c2) (0.0006935 -fitbeta8(h,a1,c1,a2,c2));

yfit = 0.0006935 -fitbeta8(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_8 covariance function','Fitted covariance')
hold off
nuggetbeta8 = 0.0006935;



%% beta9
a1 = 134.6138;
c1 = -0.0027;
a2 = 259.2714;
c2 = 0.0032;

xdata = edges';
ydata = gbeta9';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta9 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta9(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_9$', 'Interpreter', 'latex')
legend('\beta_9 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta9';
fitcovbeta9 = @(h,a1,c1,a2,c2) (0.0005396 -fitbeta9(h,a1,c1,a2,c2));

yfit = 0.0005396 -fitbeta9(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_9 covariance function','Fitted covariance')
hold off
nuggetbeta9 = 0.0005396;


%% beta10
a1 = 59.1880;
c1 = -0.0012;
a2 = 104.3672;
c2 = 0.0017;

xdata = edges';
ydata = gbeta10';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta10 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta10(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{10}$', 'Interpreter', 'latex')
legend('\beta_{10} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta10';
fitcovbeta10 = @(h,a1,c1,a2,c2) (0.0004675 -fitbeta10(h,a1,c1,a2,c2));

yfit = 0.0004675 -fitbeta10(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{10} covariance function','Fitted covariance')
hold off
nuggetbeta10 = 0.0004675;



%% beta11
a1 = 182.5439;
c1 = -0.0019;
a2 = 345.1951;
c2 = 0.0023;


xdata = edges';
ydata = gbeta11';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta11 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta11(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{11}$', 'Interpreter', 'latex')
legend('\beta_{11} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta11';
fitcovbeta11 = @(h,a1,c1,a2,c2) (0.0003845 -fitbeta11(h,a1,c1,a2,c2));

yfit = 0.0003845 -fitbeta11(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{11} covariance function','Fitted covariance')
hold off
nuggetbeta11 = 0.0003845;


%% beta12
a1 = 203.4239;%75.3961;
c1 = -0.0016;%-0.00056138;
a2 = 372.0593;%109.9840;
c2 = 0.0019;%0.00085123;


xdata = edges';
ydata = gbeta12';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta12 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta12(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{12}$', 'Interpreter', 'latex')
legend('\beta_{12} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta12';
fitcovbeta12 = @(h,a1,c1,a2,c2) (0.0003461 -fitbeta12(h,a1,c1,a2,c2));

yfit = 0.0003461 -fitbeta12(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{12} covariance function','Fitted covariance')
hold off
nuggetbeta12 = 0.0003461;




%% beta13
a1 = 117.4545;
c1 = -0.000091560;%-0.00085017;
a2 = 54.0290;%191.1814;
c2 = 0.00035336;%0.0011;


xdata = edges';
ydata = gbeta13';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta13 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta13(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{13}$', 'Interpreter', 'latex')
legend('\beta_{13} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta13';
fitcovbeta13 = @(h,a1,c1,a2,c2) (0.0002988 -fitbeta13(h,a1,c1,a2,c2));

yfit = 0.0002988 -fitbeta13(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{13} covariance function','Fitted covariance')
hold off
nuggetbeta13 = 0.0002988;





%% beta14
a1 = 93.6704;
c1 = -0.000074330;
a2 = 52.9208;
c2 = 0.00036189;


xdata = edges';
ydata = gbeta14';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta14 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta14(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{14}$', 'Interpreter', 'latex')
legend('\beta_{14} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta14';
fitcovbeta14 = @(h,a1,c1,a2,c2) (0.0003081 -fitbeta14(h,a1,c1,a2,c2));

yfit = 0.0003081 -fitbeta14(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{14} covariance function','Fitted covariance')
hold off
nuggetbeta14 = 0.0003081;





%% beta15
a1 = 126.2770;
c1 = -0.000014084;
a2 = 29.6192;
c2 = 0.00028453;


xdata = edges';
ydata = gbeta15';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta15 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta15(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{15}$', 'Interpreter', 'latex')
legend('\beta_{15} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta15';
fitcovbeta15 = @(h,a1,c1,a2,c2) (0.000275 -fitbeta15(h,a1,c1,a2,c2));

yfit = 0.000275 -fitbeta15(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{15} covariance function','Fitted covariance')
hold off
nuggetbeta15 = 0.000275;


%% beta16
a1 = 42.8983;
c1 = -0.00012339;
a2 = 55.6443;
c2 = 0.0004025;


xdata = edges';
ydata = gbeta16';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta16 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta16(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{16}$', 'Interpreter', 'latex')
legend('\beta_{16} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta16';
fitcovbeta16 = @(h,a1,c1,a2,c2) (0.0002692 -fitbeta16(h,a1,c1,a2,c2));

yfit = 0.0002692 -fitbeta16(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{16} covariance function','Fitted covariance')
hold off
nuggetbeta16 = 0.0002692;





%% beta17
a1 = 117.7087;
c1 = -0.000022107;
a2 = 28.1876;
c2 = 0.00029170;


xdata = edges';
ydata = gbeta17';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta17 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta17(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{17}$', 'Interpreter', 'latex')
legend('\beta_{17} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta17';
fitcovbeta17 = @(h,a1,c1,a2,c2) (0.0002474 -fitbeta17(h,a1,c1,a2,c2));

yfit = 0.0002474 -fitbeta17(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{17} covariance function','Fitted covariance')
hold off
nuggetbeta17 = 0.0002474;


%% beta18
a1 = 135.0416;
c1 = -0.000018408;
a2 = 26.2241;
c2 = 0.00025360;


xdata = edges';
ydata = gbeta18';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta18 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta18(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{18}$', 'Interpreter', 'latex')
legend('\beta_{18} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta18';
fitcovbeta18 = @(h,a1,c1,a2,c2) (0.0002246 -fitbeta18(h,a1,c1,a2,c2));

yfit = 0.0002246 -fitbeta18(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{18} covariance function','Fitted covariance')
hold off
nuggetbeta18 = 0.0002246;






%% beta19
a1 = 70.1403;
c1 = -0.00019985;
a2 = 31.9631;
c2 = 0.00041749;


xdata = edges';
ydata = gbeta19';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta19 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta19(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{19}$', 'Interpreter', 'latex')
legend('\beta_{19} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta19';
fitcovbeta19 = @(h,a1,c1,a2,c2) (0.0002172 -fitbeta19(h,a1,c1,a2,c2));

yfit = 0.0002172 -fitbeta19(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{19} covariance function','Fitted covariance')
hold off
nuggetbeta19 = 0.0002172;



%% beta20
a1 = 70.1510;
c1 = -0.00022249;
a2 = 31.9628;
c2 = 0.00044347;


xdata = edges';
ydata = gbeta20';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta20 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta20(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{20}$', 'Interpreter', 'latex')
legend('\beta_{20} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta20';
fitcovbeta20 = @(h,a1,c1,a2,c2) (0.0002261 -fitbeta20(h,a1,c1,a2,c2));

yfit = 0.0002261 -fitbeta20(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{20} covariance function','Fitted covariance')
hold off
nuggetbeta20 = 0.0002261;





%% beta21
a1 = 25.2853;
c1 = 0.00038791;
a2 = 45.4149;
c2 = -0.00020638;


xdata = edges';
ydata = gbeta21';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta21 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta21(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{21}$', 'Interpreter', 'latex')
legend('\beta_{21} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta21';
fitcovbeta21 = @(h,a1,c1,a2,c2) (0.0001935 -fitbeta21(h,a1,c1,a2,c2));

yfit = 0.0001935 -fitbeta21(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{21} covariance function','Fitted covariance')
hold off
nuggetbeta21 = 0.0001935;





%% beta22
a1 = 102.1439;
c1 = -0.00033842;
a2 = 63.6221;
c2 = 0.0005239;


xdata = edges';
ydata = gbeta22';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta22 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta22(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{22}$', 'Interpreter', 'latex')
legend('\beta_{22} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta22';
fitcovbeta22 = @(h,a1,c1,a2,c2) (0.0001966 -fitbeta22(h,a1,c1,a2,c2));

yfit = 0.0001966 -fitbeta22(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{22} covariance function','Fitted covariance')
hold off
nuggetbeta22 = 0.0001966;






%% beta23
a1 = 354.6484;%389.9106;%106.7431;
c1 = -0.001;%-0.0012;%-0.00049883;
a2 = 537.25;%768.3362;%77.6423;
c2 = 0.0018;%0.0024;%0.00077089;


xdata = edges';
ydata = gbeta23';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta23 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta23(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{23}$', 'Interpreter', 'latex')
legend('\beta_{23} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta23';
fitcovbeta23 = @(h,a1,c1,a2,c2) (0.0002536 -fitbeta23(h,a1,c1,a2,c2));

yfit = 0.0002536 -fitbeta23(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{23} covariance function','Fitted covariance')
hold off
nuggetbeta23 = 0.0002536;






%% beta24
a1 = 382.4297;
c1 = 0.000087988;
a2 = 27.6484;
c2 = 0.00012082;


xdata = edges';
ydata = gbeta24';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitbeta24 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+ex(h,a2,c2));  %the fitted semivariogram function
yfit = fitbeta24(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\beta_{24}$', 'Interpreter', 'latex')
legend('\beta_{24} semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfbeta24';
fitcovbeta24 = @(h,a1,c1,a2,c2) (0.0001481 -fitbeta24(h,a1,c1,a2,c2));

yfit = 0.0001481 -fitbeta24(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\beta_{24} covariance function','Fitted covariance')
hold off
nuggetbeta24 = 0.0001481;
%% Start kriging


%% Right polygon
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


%% Kriging per se
n = length(newlat);

%% Alpha1
a1 = 1003.9;
c1 = -0.1625;
a2 = 21.3262;
c2 = 0.3516;

krigged_alpha = zeros(n,1);

for p = 1:n
    krigged_alpha(p) = kriging(p,newlat, newlon, lats, lons, alpha1, fitcovalpha1, fitalpha1, nuggetalpha1, pentasph, ex, a1, c1, a2, c2);
end

krigged_alpha1 = krigged_alpha';

figure
scatter(newlon,newlat, 20, krigged_alpha1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\alpha_1$', 'Interpreter', 'latex')

save('KriggedResults_alpha1.mat','alpha1','krigged_alpha1', 'lat','lon')



%% Alpha2
a1 = 350.9542;
c1 = -0.1151;
a2 = 25.6248;
c2 = 0.2793;

krigged_alpha = zeros(n,1);

for p = 1:n
    krigged_alpha(p) = kriging(p,newlat, newlon, lats, lons, alpha2, fitcovalpha2, fitalpha2, nuggetalpha2, gaus, ex, a1, c1, a2, c2);
end

krigged_alpha2 = krigged_alpha';

figure
scatter(newlon,newlat, 20, krigged_alpha2, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\alpha_2$', 'Interpreter', 'latex')

save('KriggedResults_alpha2.mat','alpha2','krigged_alpha2', 'lat','lon')



%% Beta 1
a1 = 350.4368;
c1 = -0.1371;
a2 = 23.8057;
c2 = 0.3457;


krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta1, fitcovbeta1, fitbeta1, nuggetbeta1, gaus, ex, a1, c1, a2, c2);
end

krigged_beta1 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_1$', 'Interpreter', 'latex')

save('KriggedResults_beta1.mat','beta1','krigged_beta1', 'lat','lon')


%% Beta2
a1 = 256.6674;
c1 = -0.0025;
a2 = 4.9492;
c2 = 0.0121;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta2, fitcovbeta2, fitbeta2, nuggetbeta2, gaus, ex, a1, c1, a2, c2);
end

krigged_beta2 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta2, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_2$', 'Interpreter', 'latex')

save('KriggedResults_beta2.mat','beta2','krigged_beta2', 'lat','lon')

%% Beta3
a1 = 124.3028;
c1 = -0.0021;
a2 = 32.8914;
c2 = 0.0061;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta3, fitcovbeta3, fitbeta3, nuggetbeta3, gaus, ex, a1, c1, a2, c2);
end

krigged_beta3 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta3, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_3$', 'Interpreter', 'latex')

save('KriggedResults_beta3.mat','beta3','krigged_beta3', 'lat','lon')

%% Beta4
a1 = 98.1568;
c1 = -0.00043031;
a2 = 19.7806;
c2 = 0.0032;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta4, fitcovbeta4, fitbeta4, nuggetbeta4, gaus, ex, a1, c1, a2, c2);
end

krigged_beta4 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta4, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_4$', 'Interpreter', 'latex')

save('KriggedResults_beta4.mat','beta4','krigged_beta4', 'lat','lon')

%% Beta5 
a1 = 226.4135;
c1 = -0.0051;
a2 = 248.5329;
c2 = 0.0085;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta5, fitcovbeta5, fitbeta5, nuggetbeta5, gaus, ex, a1, c1, a2, c2);
end

krigged_beta5 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta5, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_5$', 'Interpreter', 'latex')

save('KriggedResults_beta5.mat','beta5','krigged_beta5', 'lat','lon')

%% Beta6
a1 = 180.9366;
c1 = -0.0039;
a2 = 166.5990;
c2 = 0.0058;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta6, fitcovbeta6, fitbeta6, nuggetbeta6, gaus, ex, a1, c1, a2, c2);
end

krigged_beta6 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta6, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_6$', 'Interpreter', 'latex')

save('KriggedResults_beta6.mat','beta6','krigged_beta6', 'lat','lon')


%% Beta7
a1 = -16.4444;
c1 = 0.0010;
a2 = 397.6020;
c2 = -0.00003242;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta7, fitcovbeta7, fitbeta7, nuggetbeta7, gaus, ex, a1, c1, a2, c2);
end

krigged_beta7 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta7, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_7$', 'Interpreter', 'latex')

save('KriggedResults_beta7.mat','beta7','krigged_beta7', 'lat','lon')

%% Beta8
a1 = 30.2180;%37.1070;
c1 = 0.00070866;%0.00086581;
a2 = 226.6858;%236.2554;
c2 = -0.000034576;%-0.0001985;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta8, fitcovbeta8, fitbeta8, nuggetbeta8, gaus, ex, a1, c1, a2, c2);
end

krigged_beta8 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta8, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_8$', 'Interpreter', 'latex')

save('KriggedResults_beta8.mat','beta8','krigged_beta8', 'lat','lon')


%% Beta9
a1 = 134.6138;
c1 = -0.0027;
a2 = 259.2714;
c2 = 0.0032;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta9, fitcovbeta9, fitbeta9, nuggetbeta9, gaus, sph, a1, c1, a2, c2);
end

krigged_beta9 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta9, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_9$', 'Interpreter', 'latex')

save('KriggedResults_beta9.mat','beta9','krigged_beta9', 'lat','lon')

%% Beta10
a1 = 59.1880;
c1 = -0.0012;
a2 = 104.3672;
c2 = 0.0017;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta10, fitcovbeta10, fitbeta10, nuggetbeta10, gaus, sph, a1, c1, a2, c2);
end

krigged_beta10 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta10, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{10}$', 'Interpreter', 'latex')

save('KriggedResults_beta10.mat','beta10','krigged_beta10', 'lat','lon')

%% Beta11
a1 = 182.5439;
c1 = -0.0019;
a2 = 345.1951;
c2 = 0.0023;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta11, fitcovbeta11, fitbeta11, nuggetbeta11, gaus, sph, a1, c1, a2, c2);
end

krigged_beta11 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta11, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{11}$', 'Interpreter', 'latex')

save('KriggedResults_beta11.mat','beta11','krigged_beta11', 'lat','lon')

%% Beta12
a1 = 203.4239;%75.3961;
c1 = -0.0016;%-0.00056138;
a2 = 372.0593;%109.9840;
c2 = 0.0019;%0.00085123;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta12, fitcovbeta12, fitbeta12, nuggetbeta12, gaus, sph, a1, c1, a2, c2);
end

krigged_beta12 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta12, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{12}$', 'Interpreter', 'latex')

save('KriggedResults_beta12.mat','beta12','krigged_beta12', 'lat','lon')

%% Beta13
a1 = 117.4545;
c1 = -0.000091560;%-0.00085017;
a2 = 54.0290;%191.1814;
c2 = 0.00035336;%0.0011;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta13, fitcovbeta13, fitbeta13, nuggetbeta13, gaus, sph, a1, c1, a2, c2);
end

krigged_beta13 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta13, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{13}$', 'Interpreter', 'latex')

save('KriggedResults_beta13.mat','beta13','krigged_beta13', 'lat','lon')

%% Beta14
a1 = 93.6704;
c1 = -0.000074330;
a2 = 52.9208;
c2 = 0.00036189;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta14, fitcovbeta14, fitbeta14, nuggetbeta14, gaus, sph, a1, c1, a2, c2);
end

krigged_beta14 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta14, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{14}$', 'Interpreter', 'latex')

save('KriggedResults_beta14.mat','beta14','krigged_beta14', 'lat','lon')

%% Beta15
a1 = 126.2770;
c1 = -0.000014084;
a2 = 29.6192;
c2 = 0.00028453;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta15, fitcovbeta15, fitbeta15, nuggetbeta15, gaus, sph, a1, c1, a2, c2);
end

krigged_beta15 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta15, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{15}$', 'Interpreter', 'latex')

save('KriggedResults_beta15.mat','beta15','krigged_beta15', 'lat','lon')


%% Beta16
a1 = 42.8983;
c1 = -0.00012339;
a2 = 55.6443;
c2 = 0.0004025;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta16, fitcovbeta16, fitbeta16, nuggetbeta16, gaus, sph, a1, c1, a2, c2);
end

krigged_beta16 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta16, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{16}$', 'Interpreter', 'latex')

save('KriggedResults_beta16.mat','beta16','krigged_beta16', 'lat','lon')

%% Beta17
a1 = 117.7087;
c1 = -0.000022107;
a2 = 28.1876;
c2 = 0.00029170;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta17, fitcovbeta17, fitbeta17, nuggetbeta17, gaus, ex, a1, c1, a2, c2);
end

krigged_beta17 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta17, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{17}$', 'Interpreter', 'latex')

save('KriggedResults_beta17.mat','beta17','krigged_beta17', 'lat','lon')


%% Beta18
a1 = 135.0416;
c1 = -0.000018408;
a2 = 26.2241;
c2 = 0.00025360;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta18, fitcovbeta18, fitbeta18, nuggetbeta18, gaus, ex, a1, c1, a2, c2);
end

krigged_beta18 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta18, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{18}$', 'Interpreter', 'latex')

save('KriggedResults_beta18.mat','beta18','krigged_beta18', 'lat','lon')

%% Beta19
a1 = 70.1403;
c1 = -0.00019985;
a2 = 31.9631;
c2 = 0.00041749;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta19, fitcovbeta19, fitbeta19, nuggetbeta19, gaus, ex, a1, c1, a2, c2);
end

krigged_beta19 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta19, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{19}$', 'Interpreter', 'latex')

save('KriggedResults_beta19.mat','beta19','krigged_beta19', 'lat','lon')

%% Beta20
a1 = 70.1510;
c1 = -0.00022249;
a2 = 31.9628;
c2 = 0.00044347;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta20, fitcovbeta20, fitbeta20, nuggetbeta20, gaus, ex, a1, c1, a2, c2);
end

krigged_beta20 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta20, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{20}$', 'Interpreter', 'latex')

save('KriggedResults_beta20.mat','beta20','krigged_beta20', 'lat','lon')

%% Beta21
a1 = 25.2853;
c1 = 0.00038791;
a2 = 45.4149;
c2 = -0.00020638;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta21, fitcovbeta21, fitbeta21, nuggetbeta21, gaus, ex, a1, c1, a2, c2);
end

krigged_beta21 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta21, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{21}$', 'Interpreter', 'latex')

save('KriggedResults_beta21.mat','beta21','krigged_beta21', 'lat','lon')

%% Beta22
a1 = 102.1439;
c1 = -0.00033842;
a2 = 63.6221;
c2 = 0.0005239;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta22, fitcovbeta22, fitbeta22, nuggetbeta22, gaus, ex, a1, c1, a2, c2);
end

krigged_beta22 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta22, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{22}$', 'Interpreter', 'latex')

save('KriggedResults_beta22.mat','beta22','krigged_beta22', 'lat','lon')

%% Beta23
a1 = 354.6484;%389.9106;%106.7431;
c1 = -0.001;%-0.0012;%-0.00049883;
a2 = 537.25;%768.3362;%77.6423;
c2 = 0.0018;%0.0024;%0.00077089;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta23, fitcovbeta23, fitbeta23, nuggetbeta23, gaus, ex, a1, c1, a2, c2);
end

krigged_beta23 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta23, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{23}$', 'Interpreter', 'latex')

save('KriggedResults_beta23.mat','beta23','krigged_beta23', 'lat','lon')

%% Beta24
a1 = 382.4297;
c1 = 0.000087988;
a2 = 27.6484;
c2 = 0.00012082;

krigged_beta = zeros(n,1);

for p = 1:n
    krigged_beta(p) = kriging(p,newlat, newlon, lats, lons, beta24, fitcovbeta24, fitbeta24, nuggetbeta24, gaus, ex, a1, c1, a2, c2);
end

krigged_beta24 = krigged_beta';

figure
scatter(newlon,newlat, 20, krigged_beta24, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\beta_{24}$', 'Interpreter', 'latex')

save('KriggedResults_beta24.mat','beta24','krigged_beta24', 'lat','lon')



%%
%% Alpha1
a1 = 242.6143;
a2 = 150.4824;
c1 = -0.0038;
c2 = 0.0089;
n = length(newlat);
krigged_alpha = zeros(n,1);

for p = 1:n
    krigged_alpha(p) = kriging(p,newlat, newlon, lats, lons, alpha1, fitcov1, fitalpha1, nugget1, sinehole, ex, a1, c1, a2, c2);
end

krigged_alpha1 = krigged_alpha';

figure
scatter(newlon,newlat, 20, krigged_alpha1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\alpha_1$', 'Interpreter', 'latex')

save('KriggedResults_alpha1.mat','alpha1','krigged_alpha1', 'lat','lon')


%% Alpha2
a1 = 668.6342;
a2 = 6.8191e+03;
c1 = -0.0272;
c2 = 0.2038;

n = length(newlat);
krigged_alpha = zeros(n,1);

for p = 1:n
    krigged_alpha(p) = kriging(p,newlat, newlon, lats, lons, alpha2, fitcov2, fitalpha2, nugget2, gaus, ex, a1, c1, a2, c2);
end

krigged_alpha2 = krigged_alpha';

figure
scatter(newlon,newlat, 20, krigged_alpha2, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\alpha_2$', 'Interpreter', 'latex')

save('KriggedResults_alpha2.mat','alpha2','krigged_alpha2', 'lat','lon')




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



