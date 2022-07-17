%% Kriging of the seasonality curves

% Read the parameters
 clear all
p = xlsread('NewSeasonality.xlsx', 'Hourly');

lats = p(:,2);
lons = p(:,3);
pJan = p(:,4:27);
pFeb = p(:,28:51);
pMar = p(:,52:75);
pApr = p(:,76:99);
pMay = p(:,100:123);
pJun = p(:,124:147);
pJul = p(:,148:171);
pAug = p(:,172:195);
pSep = p(:,196:219);
pOct = p(:,220:243);
pNov = p(:,244:267);
pDec = p(:,268:291);


for i=1:24
    f = char(sprintf('hour%d',i));
    v = [pJan(:,i), pFeb(:,i), pMar(:,i), pApr(:,i), pMay(:,i), pJun(:,i), pJul(:,i), pAug(:,i), pSep(:,i), pOct(:,i), pNov(:,i), pDec(:,i)];
    profile.(f) = v;
end
    
%% Make the semivariogram
for i= 1:85
    for j = 1:85
        [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
    end 
end

Ns = unique(D);
edges = 0:11:800;
[Nb, BIN] = histc(Ns, edges);

%% Build the semivariograms
for i=1:24
    f = char(sprintf('hour%d', i));
    [g.(f), cf.(f)] = calculate_CrossSemivariogram(profile.(f), edges, D);
end


%% Fit the spatial averaged semivariogram
sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
pentasph = @(x,a,c)c*(15*x./(8*a)-5/4*(x./a).^3+3/8*(x./a).^5);
gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
xdata = edges';
% for h=1:24
%     f = char(sprintf('hour%d',h));
%     for l=1:12
%         for q=1:12
%             ydata = squeeze(g.(f)(l,q,:));
%             fun1 = @(x)ssevalexp1(x,xdata,ydata);
%             x0 = [500, nanmean(ydata), 200, -5];
%             %x0 = [600,0.0004, 200, -5];
%             options = optimset('MaxFunEvals',100000);
%             bestx1 = fminsearch(fun1,x0, options);
%             a1 = bestx1(1);
%             c1 = bestx1(2);
%             a2 = bestx1(3);
%             c2 = bestx1(4);
%             yfit = pentasph(xdata,a1,c1)+gaus(xdata,a2,c2);
%                     figure;
%                     plot(xdata,ydata,'*');
%                     hold on
%                     plot(xdata,yfit,'r');
%                     xlabel('distance (km)')
%                     ylabel('Semivariogram (m/s)')
%                     title(sprintf('Cross-variogram fit for (%d,%d)',l,q))
%                     %legend('Alpha_1 semivariogram','Fitted semivariogram')
%                     hold off
%             fitprof.(f){l,q} = @(h) pentasph(h,a1,c1)+ gaus(h,a2,c2);
%         end
%     end
% end

%% Start multivariate cokringing
%save('Cross-semivariogramProfiles.mat','fitprof')
load('Cross-semivariogramProfiles.mat','fitprof')
%% Create the set of desired locations
delta = 0.1;
newlats = min(lats):0.1:max(lats);
newlons = min(lons):0.1:max(lons);
[lat, lon] = meshgrid(newlats, newlons);


newlat = [];
newlon = [];
for r=1:76
    for c=1:31
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

%% Multivariate cokring
%Calculate dist(si-s0), i=1...85
n = length(newlat); 
Dp = [];
for p=1:n
    for i=1:85
        [dist,~] = lldistkm([lats(i) lons(i)], [newlat(p) newlon(p)]);
        Dp(p,i) = dist;
    end
end

%% Hour 1
n = length(newlat); 
krigged_val = zeros(n,12);

f = char('hour1');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)', profile.(f), g.(f), fitprof.(f));
    end 
    toc
%     if p==1000
%         sprintf('Step 1000 done...')
%     elseif p==5000
%         sprintf('Step 5000 done...')
%     elseif p==10000
%         sprintf('Step 10000 done...')
%     elseif p==100000
%         sprintf('Step 100000 done...')
%     end
end


krigged_profile.(f) = krigged_val;


%% Hour 2

krigged_val = zeros(n,12);

f = char('hour2');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)', profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 3

krigged_val = zeros(n,12);

f = char('hour3');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 4

krigged_val = zeros(n,12);

f = char('hour4');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)', profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 5

krigged_val = zeros(n,12);

f = char('hour5');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)', profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 6
krigged_val = zeros(n,12);

f = char('hour6');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)', profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 7

krigged_val = zeros(n,12);

f = char('hour7');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end
    toc
   
end


krigged_profile.(f) = krigged_val;

%% Hour 8

krigged_val = zeros(n,12);

f = char('hour8');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 9

krigged_val = zeros(n,12);

f = char('hour9');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 10

krigged_val = zeros(n,12);

f = char('hour10');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 11

krigged_val = zeros(n,12);

f = char('hour11');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 12

krigged_val = zeros(n,12);

f = char('hour12');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 13

krigged_val = zeros(n,12);

f = char('hour13');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;



%% Hour 14

krigged_val = zeros(n,12);

f = char('hour14');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 15

krigged_val = zeros(n,12);

f = char('hour15');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 16

krigged_val = zeros(n,12);

f = char('hour16');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 17

krigged_val = zeros(n,12);

f = char('hour17');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 18

krigged_val = zeros(n,12);

f = char('hour18');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 19

krigged_val = zeros(n,12);

f = char('hour19');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 20

krigged_val = zeros(n,12);

f = char('hour20');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 21

krigged_val = zeros(n,12);

f = char('hour21');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Hour 22

krigged_val = zeros(n,12);

f = char('hour22');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 23

krigged_val = zeros(n,12);

f = char('hour23');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;


%% Hour 24

krigged_val = zeros(n,12);

f = char('hour24');
for p = 1:n
    tic
    for m=1:12
        krigged_val(p,m) = cokrigingVueltaTuerca(p,m,edges, newlat, newlon, lats, lons, D, Dp(p,:)',profile.(f), g.(f), fitprof.(f));
    end 
    toc
end


krigged_profile.(f) = krigged_val;

%% Re-assemble the variables
save('CokrigingResults_Seasonality.mat','newlat','newlon','profile','krigged_profile','krigged_total')




%% Make some figures




figure
scatter(newlon,newlat, 30, krigged_profile.hour1(:,12), 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for seasonality in December, hour 00')


%% Assemble the results in the right order
krigged_total = [];
for m=1:12
    for h=1:24
        f = char(sprintf('hour%d',h));
        krigged_total = [krigged_total krigged_profile.(f)(:,m)];
    end
end

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

function [g, cf] = calculate_CrossSemivariogram(z, edges, D)
%Calculate the time-space semivariogram (VerHoef&Cressie,1993)
for l=1:12
    for q=1:12
        [g(l,q,:), cf(l,q,:)] = calculate_semivariogram(z,edges,D,l,q);
    end
end

end

function [g, cf] = calculate_semivariogram(z, edges, D, h1, h2)
for k=1:length(edges)-1
    [idxrow, idxcol] = find(D>=edges(k) & D<edges(k+1));
    
    gamma = 0;
    for i=1:length(idxrow)
        gamma = gamma + (z(idxrow(i),h1)-z(idxcol(i),h1))*(z(idxrow(i),h2)-z(idxcol(i),h2));
    end
    num = length(idxrow);
    g(k) = gamma/(2*num);
    C = cov(z(idxrow,h1),z(idxcol,h2));
    cf(k) = C(1,2);
end
k = length(edges);
[idxrow, idxcol] = find(D>=edges(k-1));
gamma = 0;
for i=1:length(idxrow)
    gamma = gamma + (z(idxrow(i),h1)-z(idxcol(i),h1))*(z(idxrow(i),h2)-z(idxcol(i),h2));
end
C = cov(z(idxrow,h1),z(idxcol,h2));
cf(k) = C(1,2);
num = length(idxrow);
g(k) = gamma/(2*num);
end


