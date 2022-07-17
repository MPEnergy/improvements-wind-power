%% Kriging of the GARCH parameters

% Read the parameters
clear all
p = xlsread('GARCH.xlsx');

lats = p(:,2);
lons = p(:,3);
lambda = p(:,4);
omega = p(:,5);
gamma0 = p(:,8);
gamma1 = p(:,6);
delta = p(:,7);
skew = p(:,9);
shape = p(:,10);


%% Make the semivariogram
for i= 1:85
    for j = 1:85
        [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
    end 
end

Ns = unique(D);
edges = 0:11:800;
[Nb, BIN] = histc(Ns, edges);


%% lambda
[glambda, cflambda] = calculate_semivariogram(lambda, edges, D);
%% omega
[gomega, cfomega] = calculate_semivariogram(omega, edges, D);
%% gamma0
[ggamma0, cfgamma0] = calculate_semivariogram(gamma0, edges, D);
%% gamma1
[ggamma1, cfgamma1] = calculate_semivariogram(gamma1, edges, D);
%% delta
[gdelta, cfdelta] = calculate_semivariogram(delta, edges, D);
%% skew
[gskew, cfskew] = calculate_semivariogram(skew, edges, D);
%% shape
[gshape, cfshape] = calculate_semivariogram(shape, edges, D);


% save('Spatial_semivariogram_alphasbetas.mat','galpha1','g2','alpha1','alpha2','lats','lons','edges')

%% Fit the spatial averaged semivariogram

%% Shape
xdata = edges';
ydata = gshape';
fun1 = @(x)ssevalexp1(x,xdata,ydata);

%x0 = [500, 0.015];
 %x0 = [100, 10, 40, 1];
 x0 = [200, -0.06, 500, 0.8];
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

yfit = gaus(xdata,a1,c1)+sph(xdata,a2,c2);
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


%% lambda
a1 = 214.7325;
c1 = -0.004;
a2 = 458.9239;
c2 = 0.0053;

xdata = edges';
ydata = glambda';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitlambda = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitlambda(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\lambda$', 'Interpreter', 'latex')
legend('\lambda semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cflambda';
fitcovlambda = @(h,a1,c1,a2,c2) (0.001043 -fitlambda(h,a1,c1,a2,c2));

yfit = 0.001043 -fitlambda(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\lambda covariance function','Fitted covariance')
hold off
nuggetlambda = 0.001043;



%% omega
a1 = 472.8346;
c1 = -0.0159;
a2 = 606.6427;
c2 = 0.0154;

xdata = edges';
ydata = gomega';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitomega = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitomega(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\omega$', 'Interpreter', 'latex')
legend('\omega semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfomega';
fitcovomega = @(h,a1,c1,a2,c2) (0.00438 -fitomega(h,a1,c1,a2,c2));

yfit = 0.00438 -fitomega(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\omega covariance function','Fitted covariance')
hold off
nuggetomega = 0.00438;


%% gamma0
a1 = 328.5735;
c1 = -0.0580;
a2 = 548.5573;
c2 = 0.0616;

xdata = edges';
ydata = ggamma0';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitgamma0 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitgamma0(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\gamma_0$', 'Interpreter', 'latex')
legend('\gamma_0 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfgamma0';
fitcovgamma0 = @(h,a1,c1,a2,c2) (0.01161 -fitgamma0(h,a1,c1,a2,c2));

yfit = 0.01161 -fitgamma0(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\gamma_0 covariance function','Fitted covariance')
hold off
nuggetgamma0 = 0.01161;


%% gamma1
a1 = 333.8957;
c1 = -0.0038;
a2 = 770.3569;
c2 = 0.0054;

xdata = edges';
ydata = ggamma1';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitgamma1 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitgamma1(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\gamma_1$', 'Interpreter', 'latex')
legend('\gamma_1 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfgamma1';
fitcovgamma1 = @(h,a1,c1,a2,c2) (0.0008862 -fitgamma1(h,a1,c1,a2,c2));

yfit = 0.0008862 -fitgamma1(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\gamma_1 covariance function','Fitted covariance')
hold off
nuggetgamma1 = 0.0008862;


%% delta1
a1 = 373.5595;
c1 = -0.0256;
a2 = 565.1106;
c2 = 0.0280;

xdata = edges';
ydata = gdelta';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitdelta1 = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitdelta1(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for $\delta_1$', 'Interpreter', 'latex')
legend('\delta_1 semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfdelta';
fitcovdelta1 = @(h,a1,c1,a2,c2) (0.006838 -fitdelta1(h,a1,c1,a2,c2));

yfit = 0.006838 -fitdelta1(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('\delta_1 covariance function','Fitted covariance')
hold off
nuggetdelta1 = 0.006838;



%% skew
a1 = 294.7759;
c1 = -0.0049;
a2 = 536.0033;
c2 = 0.0073;

xdata = edges';
ydata = gskew';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitskew = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitskew(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for the skewness parameter', 'Interpreter', 'latex')
legend('Skewness semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfskew';
fitcovskew = @(h,a1,c1,a2,c2) (0.001973 -fitskew(h,a1,c1,a2,c2));

yfit = 0.001973 -fitskew(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('Skewness parameter covariance function','Fitted covariance')
hold off
nuggetskew = 0.001973;


%% shape
a1 = 171.0731;
c1 = -0.5170;
a2 = 403.6630;
c2 = 1.3667;

xdata = edges';
ydata = gshape';

sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);


fitshape = @(h,a1,c1,a2,c2) (gaus(h,a1,c1)+sph(h,a2,c2));  %the fitted semivariogram function
yfit = fitshape(xdata,a1,c1,a2,c2);
figure;

subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram')
title('Fit of spatial variation for the shape parameter', 'Interpreter', 'latex')
legend('Shape parameter semivariogram function','Fitted semivariogram')
hold off


xdata = edges';
ydata = cfshape';
fitcovshape = @(h,a1,c1,a2,c2) (0.6039 -fitshape(h,a1,c1,a2,c2));

yfit = 0.6039 -fitshape(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
legend('Shape parameter covariance function','Fitted covariance')
hold off
nuggetshape = 0.6039;

%% Start kriging


%% Right polygon
%delta = 0.01;
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


%% Lambda
a1 = 214.7325;
c1 = -0.004;
a2 = 458.9239;
c2 = 0.0053;

krigged_lambda = zeros(n,1);

for p = 1:n
    krigged_lambda(p) = kriging(p,newlat, newlon, lats, lons, lambda, fitcovlambda, fitlambda, nuggetlambda, gaus, sph, a1, c1, a2, c2);
end


figure
scatter(newlon,newlat, 20, krigged_lambda, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\lambda$', 'Interpreter', 'latex')

save('KriggedResults_lambda.mat','lambda','krigged_lambda', 'lat','lon')


%% Omega
a1 = 472.8346;
c1 = -0.0159;
a2 = 606.6427;
c2 = 0.0154;


krigged_omega = zeros(n,1);

for p = 1:n
    krigged_omega(p) = kriging(p,newlat, newlon, lats, lons, omega, fitcovomega, fitomega, nuggetomega, gaus, sph, a1, c1, a2, c2);
end


figure
scatter(newlon,newlat, 20, krigged_omega, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\omega$', 'Interpreter', 'latex')

save('KriggedResults_omega.mat','omega','krigged_omega', 'lat','lon')


%% Gamma0
a1 = 328.5735;
c1 = -0.0580;
a2 = 548.5573;
c2 = 0.0616;

krigged_gamma = zeros(n,1);

for p = 1:n
    krigged_gamma(p) = kriging(p,newlat, newlon, lats, lons, gamma0, fitcovgamma0, fitgamma0, nuggetgamma0, gaus, sph, a1, c1, a2, c2);
end

krigged_gamma0 = krigged_gamma';

figure
scatter(newlon,newlat, 20, krigged_gamma0, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\gamma_0$', 'Interpreter', 'latex')

save('KriggedResults_gamma0.mat','gamma0','krigged_gamma0', 'lat','lon')



%% Gamma1
a1 = 333.8957;
c1 = -0.0038;
a2 = 770.3569;
c2 = 0.0054;


krigged_gamma = zeros(n,1);

for p = 1:n
    krigged_gamma(p) = kriging(p,newlat, newlon, lats, lons, gamma1, fitcovgamma1, fitgamma1, nuggetgamma1, gaus, sph, a1, c1, a2, c2);
end

krigged_gamma1 = krigged_gamma';

figure
scatter(newlon,newlat, 20, krigged_gamma1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\gamma_1$', 'Interpreter', 'latex')

save('KriggedResults_gamma1.mat','gamma1','krigged_gamma1', 'lat','lon')



%% Delta1
a1 = 373.5595;
c1 = -0.0256;
a2 = 565.1106;
c2 = 0.0280;


krigged_delta = zeros(n,1);

for p = 1:n
    krigged_delta(p) = kriging(p,newlat, newlon, lats, lons, delta, fitcovdelta1, fitdelta1, nuggetdelta1, gaus, sph, a1, c1, a2, c2);
end

krigged_delta1 = krigged_delta';

figure
scatter(newlon,newlat, 20, krigged_delta1, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for $\delta_1$', 'Interpreter', 'latex')

save('KriggedResults_delta1.mat','delta','krigged_delta1', 'lat','lon')



%% Skew
a1 = 294.7759;
c1 = -0.0049;
a2 = 536.0033;
c2 = 0.0073;


krigged_skew = zeros(n,1);

for p = 1:n
    krigged_skew(p) = kriging(p,newlat, newlon, lats, lons, skew, fitcovskew, fitskew, nuggetskew, gaus, sph, a1, c1, a2, c2);
end

krigged_skew = krigged_skew';

figure
scatter(newlon,newlat, 20, krigged_skew, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for skewness parameter', 'Interpreter', 'latex')

save('KriggedResults_skew.mat','skew','krigged_skew', 'lat','lon')




%% Shape
a1 = 171.0731;
c1 = -0.5170;
a2 = 403.6630;
c2 = 1.3667;


krigged_shape = zeros(n,1);

for p = 1:n
    krigged_shape(p) = kriging(p,newlat, newlon, lats, lons, shape, fitcovshape, fitshape, nuggetshape, gaus, sph, a1, c1, a2, c2);
end

krigged_shape = krigged_shape';

figure
scatter(newlon,newlat, 20, krigged_shape, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for shape parameter', 'Interpreter', 'latex')

save('KriggedResults_shape.mat','shape','krigged_shape', 'lat','lon')


%%






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



