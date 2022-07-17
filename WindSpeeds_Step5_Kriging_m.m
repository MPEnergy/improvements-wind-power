%% Kriging of the m exponent (Weibull transform)

% Read the parameters
clear all
p = xlsread('C:\Transformations.xlsx', 'Exponent');

lats = p(:,2);
lons = p(:,3);
m = p(:,7);

%% Make the semivariogram
for i= 1:85
    for j = 1:85
        [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
    end 
end

Ns = unique(D);
edges = 0:11:800;
[Nb, BIN] = histc(Ns, edges);


%% m
z = m;
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
z = m;
gamma = 0;
for i=1:length(idxrow)
    gamma = gamma + (z(idxrow(i))-z(idxcol(i)))^2
end
C = cov(z(idxrow),z(idxcol));
cf(k) = C(1,2);
num = length(idxrow);
g(k) = gamma/(2*num);



%% Plot
figure(1);
scatter(edges, g)
xlabel('Distance(km)')
ylabel('Semivariogram')
title('Semivariogram m')



% save('Spatial_semivariogram_alphas.mat','g1','g2','alpha1','alpha2','lats','lons','edges')

%% Fit the spatial averaged semivariogram


xdata = edges';
ydata = g';
fun1 = @(x)ssevalexp1(x,xdata,ydata);

% x0 = [528, 0.0051, 528, 0.0044];
% x0 = [528, -0.003, 500, 0.0044];
x0 = [1200, 0.003, 300, 0.003];
% x0 = [50, 0.009];
options = optimset('MaxFunEvals',100000);
bestx1 = fminsearch(fun1,x0, options);

% sph = @(x,a,c0,c) c0*ones(size(x))+c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+(c0+c).*(x>a)+c0.*(x==0);
sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
ex = @(x,a,c) c.*(1-exp(-x./a));
cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
power = @(x,a,c) c.*x.^a;  %a<2, a>=0

%Check the fit
a1 = bestx1(1);
c1 = bestx1(2);
a2 = bestx1(3);
c2 = bestx1(4);

yfit = sinehole(xdata,a1,c1)+cub(xdata,a2,c2);
figure;
%subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram (m/s)')
title('Fit of spatial variation for m')
legend('Alpha_1 semivariogram','Fitted semivariogram')
hold off

sqrt(nanmean((ydata-yfit).^2))

% Fitted
a1 = 550.7927;
a2 = 502.1021;
c1 = 0.0094;
c2 = -0.0074;
fitm = @(h,a1,c1,a2,c2) (sph(h,a1,c1)+gaus(h,a2,c2));  %the fitted semivariogram function
xdata = edges;
ydata = g';

yfit = sph(xdata,a1,c1)+gaus(xdata,a2,c2);
figure;
subplot(2,1,1)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Semivariogram (m/s)')
title('Fit of spatial variation for m')
legend('m semivariogram','Fitted semivariogram')
hold off

% Fit the covariance function too
xdata = edges';
ydata = cf';

fitcovm = @(h,a1,c1,a2,c2) (0.0037-fitm(h,a1,c1,a2,c2)); 
yfit = 0.0037-fitm(xdata,a1,c1,a2,c2);
subplot(2,1,2)
plot(xdata,ydata,'*');
hold on
plot(xdata,yfit,'r');
xlabel('distance (km)')
ylabel('Covariance')
% title(sprintf('Covariance function of alpha_1'))
legend('m covariance function','Fitted covariance')
hold off
nugget = 0.0037;


sqrt(nanmean((ydata-yfit).^2))



%save('Semivariogram_m.mat','edges','g','cf', 'fitm', 'fitcovm', 'sph', 'gaus', 'a1', 'a2','c1','c2')


%% Kriging for m
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


n = length(newlat);
krigged_m = zeros(n,1);

for p = 1:n
    krigged_m(p) = kriging(p,newlat, newlon, lats, lons, m, fitcovm, fitm, nugget, sph, gaus, a1, c1, a2, c2);
end
krigged_m = krigged_m';




%%
% figure
% scatter(newlon,newlat, 20, krigged_alpha, 'filled');
% colormap(hsv) %can be removed
% xlabel('Lon')
% ylabel('Lat')
% colorbar
% title('Alpha_2')
% 
% figure
% scatter(lons,lats, 60, alpha1, 'filled');
% colormap(hsv) %can be removed
% xlabel('Lat')
% ylabel('Lon')
% colorbar
% title('Alpha_1')

%% m

figure
scatter(newlon,newlat, 20, krigged_m, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for m')

save('KriggedResults_m.mat','m','krigged_m', 'lat','lon')


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




