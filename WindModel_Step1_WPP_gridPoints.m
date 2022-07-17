%% Step1: Build the wind power prediction model in-sample with the nearby grid points
WSin = load('GridData.mat');
WSoos = load('GridData_oos.mat');
datesin = WSin.dat;
valuesin = WSin.val;
datesoos = WSoos.dat;
valuesoos = WSoos.val;
lats = WSin.lats;
lons = WSin.lons;
%% Reading in-sample data
% Point 1: (35.5,-118)
idx = find(lats==35.5 & lons==-118);
ws1 = valuesin(:,idx);
%Point 2: (35, -118)
idx = find(lats==35 & lons==-118);
ws2 = valuesin(:,idx);
%Point 3: (35, -118.5)
idx = find(lats==35 & lons==-118.5);
ws3 = valuesin(:,idx);

%Point 4: (34,-116.5)
idx = find(lats==34 & lons==-116.5);
ws4 = valuesin(:,idx);

%Point 5: (33,-115.5)
idx = find(lats==33 & lons==-115.5);
ws5 = valuesin(:,idx);
%Point 6: (33, -116.5)
idx = find(lats==33 & lons==-116.5);
ws6 = valuesin(:,idx);

%% Reading out-of-sample data
% Point 1: (35.5,-118)
idx = find(lats==35.5 & lons==-118);
ws1oos = valuesoos(:,idx);
%Point 2: (35, -118)
idx = find(lats==35 & lons==-118);
ws2oos = valuesoos(:,idx);
%Point 3: (35, -118.5)
idx = find(lats==35 & lons==-118.5);
ws3oos = valuesoos(:,idx);

%Point 4: (34,-116.5)
idx = find(lats==34 & lons==-116.5);
ws4oos = valuesoos(:,idx);

%Point 5: (33,-115.5)
idx = find(lats==33 & lons==-115.5);
ws5oos = valuesoos(:,idx);
%Point 6: (33, -116.5)
idx = find(lats==33 & lons==-116.5);
ws6oos = valuesoos(:,idx);




%% Wind farm specific information:
% Locations of wind farms
lat1 = 35.2864;
lon1 = -118.1948;
cap1 = [370 370 370 370 370 370];

lat2 = [35.0634 35.0619 35.0069];
lon2 = [-118.3742 -118.2931 -118.2246];
cap2 = [2094 2094 2094 2287 2418 2490];

lat3 = 34.9095;
lon3 = -118.4410;
cap3 = [333 333 333 333 333 333];

lat4 = [33.9138];
lon4 = [-116.5871];
cap4 = [606 608 608 609 609 648];

lat5 = [32.7478];
lon5 = [-116.0717];
cap5 = [264 264 264 264 264 264];

lat6 = [32.7471];
lon6 = [-116.2829];
cap6 = [52 52 52 183 183 183];


%% Wind Power generation in SP15
WPP = load('WPP.mat');


cap = [cap1' cap2' cap3' cap4' cap5' cap6'];
%cap = cap1+cap2+cap3+cap4+cap5+cap6;
capTotal = [3697 3698 3700 4026 4157 4268]; %some small wind mills should be added
dat = [datenum('1-Jan-2015'), datenum('1-Jan-2016'), datenum('1-Jan-2017'),datenum('1-Jan-2018'), datenum('1-Jan-2019'), datenum('2-Mar-2020')];
cap = interp1(dat,cap,WPP.dat,'linear');
CF = WPP.val./cap;
CF(CF<0) = 0;


% In-sample market data
[datesin, i1, i2] = intersect(round(datesin.*24)./24, round(WPP.dat.*24)./24);
prodin = WPP.val(i2);
prodin(prodin<0) = 0;
capin = cap(i2,:);
CFin = CF(i2);
ws1 = ws1(i1);
ws2 = ws2(i1);
ws3 = ws3(i1);
ws4 = ws4(i1);
ws5 = ws5(i1);
ws6 = ws6(i1);
capTotalin = sum(capin,2);

% Out-of-sample market data
[datesoos, i1, i2] = intersect(round(datesoos.*24)./24, round(WPP.dat.*24)./24);
prodoos = WPP.val(i2);
prodoos(prodoos<0) = 0;
capoos = cap(i2,:);
CFoos = CF(i2);
ws1oos = ws1oos(i1);
ws2oos = ws2oos(i1);
ws3oos = ws3oos(i1);
ws4oos = ws4oos(i1);
ws5oos = ws5oos(i1);
ws6oos = ws6oos(i1);
capTotaloos = sum(capoos,2);


%% Build the wind model per se
% 1. Wind speeds for each station
v = [ws1 ws2 ws3 ws4 ws5 ws6];

%% 2. Rotor dimensions: A
diameter = [90 80 92.5 61.5 108 107];
A = (diameter./2).^2.*pi;

%% 3. Cut-in speeds
cutin = [3.5 3 3 2.5 3 3];

%% 4. Cut-out speeds
cutout = [25 25 24 25 25 23];

%% 5. Rated speeds
rated = [11 13 12.5 13 11.5 14];
%% 6. Fourier terms for seasonality of air density
f = [];
dates = datesin;
for st=1:6
    numCycles = 2;
    vv = round(dates.*24)./24-datenum(year(dates),1,1);
    cycles = repmat(struct('dat',[],'name',[]),1,2*numCycles);
    for ii = 1:numCycles
        tc = cos(2*pi*vv/365.25*ii);
        ts = sin(2*pi*vv/365.25*ii);
        cycles(2*(ii-1)+1).dat = tc;
        cycles(2*(ii-1)+1).name = ['cos' num2str(ii)];
        cycles(2*ii).dat = ts;
        cycles(2*ii).name = ['sin' num2str(ii)];
    end
    
    
    numCycles2 = 0;
    vv = hour(dates);
    cyclesH = repmat(struct('dat',[],'name',[]),1,2*numCycles2);
    for ii = 1:numCycles2
        tc = cos(2*pi*vv/24*ii);
        ts = sin(2*pi*vv/24*ii);
        cyclesH(2*(ii-1)+1).dat = tc;
        cyclesH(2*(ii-1)+1).name = ['cosH' num2str(ii)];
        cyclesH(2*ii).dat = ts;
        cyclesH(2*ii).name = ['sinH' num2str(ii)];
    end
   f = [f ones(size([cycles.dat],1), 1) cycles.dat cyclesH.dat];
end

%% Start the calibration process
options = optimset('MaxFunEvals', 5000, 'MaxIter', 10000, 'TolFun', 1e-6, 'TolX', 1e-6, 'TolCon', 1e-6, 'AlwaysHonorConstraints', 'bounds', 'Display', 'iter');%, 'Algorithm','sqp');


numFarms = 6;
estimatedProd = nan(size(dates));

% tart taking 1000 random weights and choose the ones giving the best error
w = 1334.7*rand(1000, numFarms)+71.6;
errorWPP = local_calculateWPPerror(w, v,cutin, cutout, rated,  capin, repmat(capTotalin',1000,1),  A, f, repmat(prodin',1000,1), 1);

[~,I] = sort(errorWPP);
w0 = w(I(1),:);
x0 = w0';
 
 % The bounds are decided as approx #turbines in wind park*C_f Betz' limit
 % numbers
LB=[repmat(71.6,size(x0,1),1)];       
UB=[repmat(1406.3,size(x0,1),1)];
% the opitmization process 
[paramsFit,fv] = fminsearchbnd(@(x) local_calculateWPPerror(x, v, cutin, cutout, rated,  capin, capTotalin, A, f, prodin', 0), ...
        x0, LB,UB, options);
 
w = paramsFit';
estimatedProd = local_calculateWPP(w, v,cutin, cutout, rated,  capin, A, f);
estimatedProd(estimatedProd>capTotalin) = capTotalin(estimatedProd>capTotalin);
 
figure;
hold on;
tSeries(dates,estimatedProd).plot('r')
tSeries(dates, prodin).plot
xlabel('Date')
ylabel('MWh')
title('In-sample Wind Power Prediction with grid points')


MAE = mean(abs(estimatedProd-prodin))
%MAPE = mean(abs(estimatedProd(prodin~=0)-prodin(prodin~=0))./abs(prodin(prodin~=0)))
MAPE = 100*mean(abs(estimatedProd(prodin~=0)-prodin(prodin~=0)))/mean(prodin(prodin~=0))
RMSE = sqrt(mean((estimatedProd-prodin).^2))

%% Compare to the naive model
% fcastedProd(1) = prodin(1);
% for i=2:length(prodin)
%     fcastedProd(i) = prodin(i-1);
% end
% fcastedProd = fcastedProd';
% 
% MAE = mean(abs(fcastedProd-prodin))
% MAPE =100* mean(abs(fcastedProd(prodin~=0)-prodin(prodin~=0))./mean(prodin(prodin~=0)))
% RMSE = sqrt(mean((fcastedProd-prodin).^2))






%% Check the out-of-sample
voos = [ws1oos ws2oos ws3oos ws4oos ws5oos ws6oos];
estimatedProdoos = local_calculateWPP(w, voos,cutin, cutout, rated,  capoos, A, f);

figure;
hold on;
tSeries(datesoos,estimatedProdoos).plot('r')
tSeries(datesoos, prodoos).plot
xlabel('Date')
ylabel('MWh')
title('Out-of-sample Wind Power Prediction with grid points')


MAE = mean(abs(estimatedProdoos-prodoos))
MAPE = 100*mean(abs(estimatedProdoos(prodoos~=0)-prodoos(prodoos~=0)))/mean(prodoos(prodoos~=0))
RMSE = sqrt(mean((estimatedProdoos-prodoos).^2))


%% Code: target function
function errorWPP = local_calculateWPPerror(params,v, cutin, cutout, rated, cap, capTotal, A, f, prod, multiple)
% Function to estimate the weighted power curves of the separate wind farms
estimPower = 0;
numFarms = length(A);
if multiple
    w = params;
else
    w = params;
end


for i=1:numFarms
    if multiple
        
        estimPower = estimPower+((((w(i)*(1/2)*A(i)*1.06).*v(:,i)'.^3))).*(v(:,i)'>=cutin(i) & v(:,i)'<rated(i))+...
            0.593.*(cap(:,i)'.*10^3).*(v(:,i)'>=rated(i) & v(:,i)'<cutout(i));
    else
        
        estimPower = estimPower+(((w(i)*(1/2)*A(i)*1.06).*v(:,i).^3)).*(v(:,i)>=cutin(i) & v(:,i)<rated(i))+...
            0.593.*(cap(:,i).*10^3).*(v(:,i)>=rated(i) & v(:,i)<cutout(i));
    end
    
    
end
 estimPower= estimPower./10^6; 
if multiple
    errorWPP = mean((prod-estimPower).^2./capTotal.^2,2);
else
    errorWPP = double(mean((prod'-estimPower).^2./capTotal.^2));
end

end


function estimPower = local_calculateWPP(params,v, cutin, cutout, rated, cap, A, f)
% Function to estimate the weighted power curves of the separate wind farms
estimPower = 0;
numFarms = length(A);
w = params;

for i=1:numFarms
  
 estimPower = estimPower+(((w(i)*(1/2)*A(i)*1.06).*v(:,i).^3)).*(v(:,i)>=cutin(i) & v(:,i)<rated(i))+...
        0.593.*(cap(:,i).*10^3).*(v(:,i)>=rated(i) & v(:,i)<cutout(i));

end

estimPower = estimPower./10^6;

end
