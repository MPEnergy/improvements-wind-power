%% Step3: Calibrate the model in-sample for the grid points, 
%do the same for the wind farm locations. Then use the same parameters to derive out-of-sample predictions
clear all
% Read the data
WSin = load('GridData.mat');
WSoos = load('GridData_oos.mat');
datesin = WSin.dat;
valuesin = WSin.val;
datesoos = WSoos.dat;
valuesoos = WSoos.val;
lats = WSin.lats;
lons = WSin.lons;

%Read the krigged data
kriggedin2015 = load('Results_Kriging_WindSpeeds2015.mat');
kriggedin2016 = load('Results_Kriging_WindSpeeds2016.mat');
kriggedin2017 = load('Results_Kriging_WindSpeeds2017.mat');
kriggedin2018 = load('Results_Kriging_WindSpeeds2018.mat');
kriggedin2019 = load('Results_Kriging_WindSpeeds2019.mat');
kriggedoos = load('Results_Kriging_OOS_WindSpeeds.mat','datesoos', 'krigged_result');
%kriggedpureoos = load('Results_Kriging_PureOOS_WindSpeeds.mat','datesoos', 'krigged_result');

load('NewCoords.mat','newlat','newlon')

%% Wind farm specific data
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

% lat4 = [33.9523 33.9156 33.9036];%[33.9138];
% lon4 = [-116.6589 -116.7172 -116.5694];%[-116.5871];
% cap4 = [56 56 56 56 56 56; 85 85 85 85 85 85; 431 432 432 432 432 471];%[606 608 608 609 609 648];
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

%% Wind Speed in Grid Point:In-sample market data

%Find the right ws
% Point 1: (35.5,-118)
idx = find(newlat==round(35.5*100/100,2) & newlon==round(-118*100/100,2)); 
ws1 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
%Point 2: (35, -118)
idx = find(newlat==round(35*100/100,2) & newlon==round(-118*100/100,2)); 
ws2 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
%Point 3: (35, -118.5)
idx = find(newlat==round(35*100/100,2) & newlon==round(-118.5*100/100,2)); 
ws3 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
%Point 4: (34,-116.5)
idx = find(newlat==round(34*100/100,2) & newlon==round(-116.5*100/100,2)); 
ws4 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
%Point 5: (33,-115.5)
idx = find(newlat==round(33*100/100,2) & newlon==round(-115.5*100/100,2)); 
ws5 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
%Point 6: (33, -116.5)
idx = find(newlat==round(33*100/100,2) & newlon==round(-116.5*100/100,2)); 
ws6 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

[datesin, i1, i2] = intersect(round(datesin.*24)./24, round(WPP.dat.*24)./24);
prodin = WPP.val(i2);
prodin(prodin<0) = 0;
capin = cap(i2,:);
ws1 = ws1(i1);
ws2 = ws2(i1);
ws3 = ws3(i1);
ws4 = ws4(i1);
ws5 = ws5(i1);
ws6 = ws6(i1);

% Out-of-sample in grid points

% Point 1: (35.5,-118)
idx = find(newlat==round(35.5*100/100,2) & newlon==round(-118*100/100,2)); 
ws1oos = kriggedoos.krigged_result(idx,:);
%Point 2: (35, -118)
idx = find(newlat==round(35*100/100,2) & newlon==round(-118*100/100,2)); 
ws2oos = kriggedoos.krigged_result(idx,:);
%Point 3: (35, -118.5)
idx = find(newlat==round(35*100/100,2) & newlon==round(-118.5*100/100,2)); 
ws3oos = kriggedoos.krigged_result(idx,:);
%Point 4: (34,-116.5)
idx = find(newlat==round(34*100/100,2) & newlon==round(-116.5*100/100,2)); 
ws4oos = kriggedoos.krigged_result(idx,:);
%Point 5: (33,-115.5)
idx = find(newlat==round(33*100/100,2) & newlon==round(-115.5*100/100,2)); 
ws5oos = kriggedoos.krigged_result(idx,:);
%Point 6: (33, -116.5)
idx = find(newlat==round(33*100/100,2) & newlon==round(-116.5*100/100,2)); 
ws6oos = kriggedoos.krigged_result(idx,:);

[datesoos, i1oos, i2oos] = intersect(round(datesoos.*24)./24, round(WPP.dat.*24)./24);
prodoos = WPP.val(i2oos);
prodoos(prodoos<0) = 0;
capoos = cap(i2oos,:);
ws1oos = ws1oos(i1oos);
ws2oos = ws2oos(i1oos);
ws3oos = ws3oos(i1oos);
ws4oos = ws4oos(i1oos);
ws5oos = ws5oos(i1oos);
ws6oos = ws6oos(i1oos);

%% Wind Speed at Wind Farm: In-sample market data

%Find the right ws
idx = find(newlat==round(lat1*100/100,2) & newlon==round(lon1*100/100,2)); 
wf1 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat2(1)*100/100,2) & newlon==round(lon2(1)*100/100,2)); 
wf21 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat2(2)*100/100,2) & newlon==round(lon2(2)*100/100,2)); 
wf22 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat2(3)*100/100,2) & newlon==round(lon2(3)*100/100,2)); 
wf23 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat3*100/100,2) & newlon==round(lon3*100/100,2)); 
wf3 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat4*100/100,2) & newlon==round(lon4*100/100,2)); 
wf4 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
% idx = find(newlat==round(lat4(1)*100/100,2) & newlon==round(lon4(1)*100/100,2)); 
% wf41 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
%     kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
% 
% idx = find(newlat==round(lat4(2)*100/100,2) & newlon==round(lon4(2)*100/100,2)); 
% wf42 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
%     kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';
% 
% idx = find(newlat==round(lat4(3)*100/100,2) & newlon==round(lon4(3)*100/100,2)); 
% wf43 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
%     kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat5*100/100,2) & newlon==round(lon5*100/100,2)); 
wf5 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

idx = find(newlat==round(lat6*100/100,2) & newlon==round(lon6*100/100,2)); 
wf6 = [kriggedin2015.krigged_result2015(idx,:), kriggedin2016.krigged_result2016(idx,:), kriggedin2017.krigged_result2017(idx,:),...
    kriggedin2018.krigged_result2018(idx,:), kriggedin2019.krigged_result2019(idx,:)]';

% [datesin, i1, i2] = intersect(round(datesin.*24)./24, round(WPP.dat.*24)./24);
% prodin = WPP.val(i2);
% prodin(prodin<0) = 0;
% capin = cap(i2,:);
wf1 = wf1(i1);
wf21 = wf21(i1);
wf22 = wf22(i1);
wf23 = wf23(i1);
wf3 = wf3(i1);
wf4 = wf4(i1);
% wf41 = wf41(i1);
% wf42 = wf42(i1);
% wf43 = wf43(i1);
wf5 = wf5(i1);
wf6 = wf6(i1);
capTotalin = sum(capin,2);

% Out-of-sample market data

idx = find(newlat==round(lat1*100/100,2) & newlon==round(lon1*100/100,2)); 
wf1oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat2(1)*100/100,2) & newlon==round(lon2(1)*100/100,2)); 
wf21oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat2(2)*100/100,2) & newlon==round(lon2(2)*100/100,2)); 
wf22oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat2(3)*100/100,2) & newlon==round(lon2(3)*100/100,2)); 
wf23oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat3*100/100,2) & newlon==round(lon3*100/100,2)); 
wf3oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat4*100/100,2) & newlon==round(lon4*100/100,2)); 
wf4oos = kriggedoos.krigged_result(idx,:);
% idx = find(newlat==round(lat4(1)*100/100,2) & newlon==round(lon4(1)*100/100,2)); 
% wf41oos = kriggedoos.krigged_result(idx,:);
% 
% idx = find(newlat==round(lat4(2)*100/100,2) & newlon==round(lon4(2)*100/100,2)); 
% wf42oos = kriggedoos.krigged_result(idx,:);
% 
% idx = find(newlat==round(lat4(3)*100/100,2) & newlon==round(lon4(3)*100/100,2)); 
% wf43oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat5*100/100,2) & newlon==round(lon5*100/100,2)); 
wf5oos = kriggedoos.krigged_result(idx,:);

idx = find(newlat==round(lat6*100/100,2) & newlon==round(lon6*100/100,2)); 
wf6oos = kriggedoos.krigged_result(idx,:);

% [datesoos, i1, i2] = intersect(round(datesoos.*24)./24, round(WPP.dat.*24)./24);
% prodoos = WPP.val(i2);
% prodoos(prodoos<0) = 0;
% capoos = cap(i2,:);
wf1oos = wf1oos(i1oos);
wf21oos = wf21oos(i1oos);
wf22oos = wf22oos(i1oos);
wf23oos = wf23oos(i1oos);
wf3oos = wf3oos(i1oos);
wf4oos = wf4oos(i1oos);
% wf41oos = wf41oos(i1oos);
% wf42oos = wf42oos(i1oos);
% wf43oos = wf43oos(i1oos);
wf5oos = wf5oos(i1oos);
wf6oos = wf6oos(i1oos);
capTotaloos = sum(capoos,2);

%% Build the wind model per se

%% 1. Calibrate in-sample & measure out-of-sample in grid points
% Rotor dimensions: A
diameter = [90 80 92.5 61.5 108 107];
A = (diameter./2).^2.*pi;
% Cut-in speeds
%cutin = [3.5 3 3 2.5 3 3];
cutin = [1.5 1 1.5 1 1 1];

% Cut-out speeds
cutout = [25 24.33 24 25 25 23];
% Rated speeds
%rated = [11 13 12.5 13 11.5 14];
%%rated = [11 13 12.5 13 10 12.5];
rated = [11 8.5 12.5 11 8 8.5];%rated = [11 8.5 8.5 8.5 12.5 10 11 9 8 8.5];
%Actual in-sample Wind speeds in the grid points
v = [ws1 ws2 ws3 ws4 ws5 ws6];
% Forecasted in-sample Wind speeds in the grid points
voos = [ws1oos' ws2oos' ws3oos' ws4oos' ws5oos' ws6oos'];
%capin2 = [capin(:,1) capin(:,2) capin(:,3) capin(:,4)+capin(:,5)+capin(:,6) capin(:,7) capin(:,8)];
%Load the randomness factor
load('RecreateRandomness.mat')
% rng('default')
% s = rng;
[w_gridPointsActual, estimatedProd_gridPoints_insample] = util_CalibrateWPP(datesin, v, cutin, cutout, rated, A, capin, capTotalin, prodin, s);
estimatedProd_gridPoints_oos = local_calculateWPP(w_gridPointsActual, voos, cutin, cutout, rated,  capoos, A, []);

sprintf('Grid Points in-sample errors:')
MAE = mean(abs(estimatedProd_gridPoints_insample-prodin))
MAPE = 100*mean(abs(estimatedProd_gridPoints_insample(prodin~=0)-prodin(prodin~=0)))/mean(prodin(prodin~=0))
RMSE = sqrt(mean((estimatedProd_gridPoints_insample-prodin).^2))

sprintf('Grid Points out-of-sample errors:')
MAE = mean(abs(estimatedProd_gridPoints_oos-prodoos))
MAPE = 100*mean(abs(estimatedProd_gridPoints_oos(prodoos~=0)-prodoos(prodoos~=0)))/mean(prodoos(prodoos~=0))
RMSE = sqrt(mean((estimatedProd_gridPoints_oos-prodoos).^2))



figure;
hold on;
tSeries(datesin,estimatedProd_gridPoints_insample).plot('-*k')
tSeries(datesoos,estimatedProd_gridPoints_oos).plot('-*r')
tSeries(datesin, prodin).plot('k')
tSeries(datesoos, prodoos).plot('r')
xlabel('Date')
ylabel('MWh')
title('Wind Power Prediction with gridded wind speeds')



%% 3. Calibrate in-sample & measure out-of-sample with Wind speeds at wind farm locations
% Rotor dimensions: A
%diameter = [90 90 80 80 92.5 61.5 108 107];
%diameter = [90 90 80 80 92.5 15.3 47 45 108 107];
diameter = [90 90 80 80 92.5 61.5 108 107];%[90 90 80 80 92.5 15.3 47 45 108 107];%diameter = [90 90 80 80 92.5 61.5 108 107];%diameter = [90 80 92.5 61.5 108 107];
A = (diameter./2).^2.*pi;
% Cut-in speeds
%cutin = [3.5 3 3 3 3 2.5 3 3];
%cutin = [3.5 2.5 2.5 2.5 3 2.5 3 3];
%cutin = [3.5 2.5 2.5 2.5 3 2.5 2.5 2.5 3 3];
% cutin = [1.5 1 1 1 1.5 1 1 1 1 1];
cutin = [1.5 1 1 1 1.5 1 1 1];%[1.5 1 1 1 1.5 1 1 1 1 1];%cutin = [1.5 1 1 1 1.5 1 1 1];%cutin = [1.5 1 1.5 1 1 1];
% Cut-out speeds
%cutout = [25 25 25 23 24 25 25 23];
%cutout = [25 25 25 23 24 25 25 25 25 23];
cutout =  [25 25 25 23 24 25 25 23];%[25 25 25 23 24 25 25 25 25 23];%cutout = [25 25 25 23 24 25 25 23];%cutout = [25 25 24 25 25 23];
% Rated speeds
%rated = [11 15 13 12.5 12.5 13 11.5 14];
%rated = [11 15 12 10 12.5 13 11.5 14];
%rated = [11 15 12 10 12.5 16 11.3 13 11.5 14];
%rated = [11 8.5 8.5 8.5 12.5 10 11 9 8 8.5];
rated = [11 8.5 8.5 8.5 12.5 11 8 8.5];%[11 8.5 8.5 8.5 12.5 10 11 9 8 8.5];%rated = [11 8.5 8.5 8.5 12.5 11 8 8.5];%rated = [11 8.5 12.5 11 8 8.5];
% Capacities
% capin = [capin(:,1) capin(:,2) capin(:,3) capin(:,4) capin(:,5) capin(:,6)];
capin = [capin(:,1) capin(:,2)./3 capin(:,2)./3 capin(:,3)./3 capin(:,3) capin(:,4) capin(:,5) capin(:,6)];
%capin = [capin(:,1) capin(:,2)./3 capin(:,2)./3 capin(:,2)./3 capin(:,3) capin(:,4) capin(:,5) capin(:,6) capin(:,7) capin(:,8)];

% vWF = [wf1 wf21 wf3 wf4 wf5 wf6];
vWF = [wf1 wf21 wf22 wf23 wf3 wf4 wf5 wf6];
% vWF = [wf1 wf21 wf22 wf23 wf3 wf41 wf42 wf43 wf5 wf6];
% vWF = [wf1 wf21 wf22 wf23 wf3 wf4 wf5 wf6];
% Out-of-sample
% capoos = [capoos(:,1) capoos(:,2) capoos(:,3) capoos(:,4) capoos(:,5) capoos(:,6)];
capoos = [capoos(:,1) capoos(:,2)./3 capoos(:,2)./3 capoos(:,2)./3  capoos(:,3) capoos(:,4) capoos(:,5) capoos(:,6)];
%capoos = [capoos(:,1) capoos(:,2)./3 capoos(:,2)./3 capoos(:,2)./3 capoos(:,3) capoos(:,4) capoos(:,5) capoos(:,6) capoos(:,7) capoos(:,8)];
% vWFoos = [wf1oos' wf21oos' wf3oos' wf4oos' wf5oos' wf6oos'];
vWFoos = [wf1oos' wf21oos' wf22oos' wf23oos' wf3oos' wf4oos' wf5oos' wf6oos'];
% vWFoos = [wf1oos' wf21oos' wf22oos' wf23oos' wf3oos' wf41oos' wf42oos' wf43oos' wf5oos' wf6oos'];

[w_windFarm, estimatedProd_WindFarm_insample] = util_CalibrateWPP(datesin, vWF, cutin, cutout, rated, A, capin, capTotalin, prodin, s);
estimatedProd_WindFarm_oos = local_calculateWPP(w_windFarm, vWFoos,cutin, cutout, rated,  capoos, A, []);

sprintf('In-sample at wind farm locations:')
MAE = mean(abs(estimatedProd_WindFarm_insample-prodin))
MAPE = 100*mean(abs(estimatedProd_WindFarm_insample(prodin~=0)-prodin(prodin~=0)))/mean(prodin(prodin~=0))
RMSE = sqrt(mean((estimatedProd_WindFarm_insample-prodin).^2))


sprintf('Out-of-sample at wind farm locations:')

MAE = mean(abs(estimatedProd_WindFarm_oos-prodoos))
MAPE = 100*mean(abs(estimatedProd_WindFarm_oos(prodoos~=0)-prodoos(prodoos~=0)))/mean(prodoos(prodoos~=0))
RMSE = sqrt(mean((estimatedProd_WindFarm_oos-prodoos).^2))



figure;
hold on;
tSeries(datesin,estimatedProd_WindFarm_insample).plot('-*k')
tSeries(datesoos,estimatedProd_WindFarm_oos).plot('-*r')
tSeries(datesin, prodin).plot('k')
tSeries(datesoos, prodoos).plot('r')
xlabel('Date')
ylabel('MWh')
title('Wind Power Prediction with onsite wind speeds')


%% Calculate confidence intervals
gridPointsp = estimatedProd_gridPoints_insample+1.96*std(estimatedProd_gridPoints_insample);
gridPointsm = estimatedProd_gridPoints_insample-1.96*std(estimatedProd_gridPoints_insample);

idxGridPoints = find(prodin<gridPointsm | prodin>gridPointsp);

sprintf('In-sample Grid Points: %d%% outside the CI',length(idxGridPoints)*100/length(prodin))


X=[datesin;flipud(datesin)];               
Y=[gridPointsm;flipud(gridPointsp)]; 
figure; hold on;
fill(X,Y,[0.7    0.7    0.7]);
scatter(datesin, prodin,'r')
tSeries(datesin, gridPointsm).plot('k');
tSeries(datesin, gridPointsp).plot('k');
title('In-sample model based on grided data: 95-5% PI')








gridPointsoosp = estimatedProd_gridPoints_oos+1.96*std(estimatedProd_gridPoints_oos);
gridPointsoosm = estimatedProd_gridPoints_oos-1.96*std(estimatedProd_gridPoints_oos);

idxGridPointsoos = find(prodoos<gridPointsoosm | prodoos>gridPointsoosp);

sprintf('Out-of-sample Grid Points: %d%% outside the CI',length(idxGridPointsoos)*100/length(prodoos))

X=[datesoos;flipud(datesoos)];               
Y=[gridPointsoosm;flipud(gridPointsoosp)]; 
figure; hold on;
fill(X,Y,[0.7    0.7    0.7]);
scatter(datesoos, prodoos,'r')
tSeries(datesoos, gridPointsoosm).plot('k');
tSeries(datesoos, gridPointsoosp).plot('k');
title('Out-of-sample model based on grided data: 95-5% PI')





windFarmp = estimatedProd_WindFarm_insample+1.96*std(estimatedProd_WindFarm_insample);
windFarmm = estimatedProd_WindFarm_insample-1.96*std(estimatedProd_WindFarm_insample);

idxWindFarm = find(prodin<windFarmm | prodin>windFarmp);

sprintf('In-sample Wind Farm Points: %d%% outside the CI',length(idxWindFarm)*100/length(prodin))

X=[datesin;flipud(datesin)];               
Y=[windFarmm;flipud(windFarmp)]; 
figure; hold on;
fill(X,Y,[0.7    0.7    0.7]);
scatter(datesin, prodin,'r')
tSeries(datesin, windFarmm).plot('k');
tSeries(datesin, windFarmp).plot('k');
title('In-sample model based on kriged data: 95-5% PI')




windFarmoosp = estimatedProd_WindFarm_oos+1.96*std(estimatedProd_WindFarm_oos);
windFarmoosm = estimatedProd_WindFarm_oos-1.96*std(estimatedProd_WindFarm_oos);

idxWindFarmoos = find(prodoos<windFarmoosm | prodoos>windFarmoosp);

sprintf('Out-of-sample Wind Farm Points: %d%% outside the CI',length(idxWindFarmoos)*100/length(prodoos))

X=[datesoos;flipud(datesoos)];               
Y=[windFarmoosm;flipud(windFarmoosp)]; 
figure; hold on;
fill(X,Y,[0.7    0.7    0.7]);
scatter(datesoos, prodoos,'r')
tSeries(datesoos, windFarmoosm).plot('k');
tSeries(datesoos, windFarmoosp).plot('k');
title('Out-of-sample model based on kriged data: 95-5% PI')





























function [w, estimatedProd] = local_CalibrateWPP(dates, v, cutin, cutout, rated, A, cap, capTotal, prod)

options = optimset('MaxFunEvals', 5000, 'MaxIter', 10000, 'TolFun', 1e-6, 'TolX', 1e-6, 'TolCon', 1e-6, 'AlwaysHonorConstraints', 'bounds', 'Display', 'iter');%, 'Algorithm','sqp');


numFarms = size(v,2);
estimatedProd = nan(size(dates));
f = [];
% Start taking 1000 random weights and choose the ones giving the best error
w = 1334.7*rand(1000, numFarms)+71.6;
errorWPP = local_calculateWPPerror(w, v,cutin, cutout, rated,  cap, repmat(capTotal',1000,1),  A, f, repmat(prod',1000,1), 1);

[~,I] = sort(errorWPP);
w0 = w(I(1),:);
x0 = w0';

% The bounds are decided as approx #turbines in wind park*C_f Betz' limit
% numbers
LB=[repmat(71.6,size(x0,1),1)];
UB=[repmat(1406.3,size(x0,1),1)];
% the opitmization process
[paramsFit,fv] = fminsearchbnd(@(x) local_calculateWPPerror(x, v, cutin, cutout, rated,  cap, capTotal, A, f, prod', 0), ...
    x0, LB,UB, options);

w = paramsFit';
estimatedProd = local_calculateWPP(w, v,cutin, cutout, rated,  cap, A, f);
estimatedProd(estimatedProd>capTotal) = capTotal(estimatedProd>capTotal);


end





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
        
        estimPower = estimPower+((((w(:,i)*(1/2)*A(i)*1.225).*v(:,i)'.^3))).*(v(:,i)'>=cutin(i) & v(:,i)'<rated(i))+...
            0.593.*(cap(:,i)'.*10^3).*(v(:,i)'>=rated(i) & v(:,i)'<cutout(i));
    else
        
        estimPower = estimPower+(((w(i)*(1/2)*A(i)*1.225).*v(:,i).^3)).*(v(:,i)>=cutin(i) & v(:,i)<rated(i))+...
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
  
 estimPower = estimPower+(((w(i)*(1/2)*A(i)*1.225).*v(:,i).^3)).*(v(:,i)>=cutin(i) & v(:,i)<rated(i))+...
        0.593.*(cap(:,i).*10^3).*(v(:,i)>=rated(i) & v(:,i)<cutout(i));

end

estimPower = estimPower./10^6;

end









