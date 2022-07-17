%% Step2: Build wind speed predictions 1h-6h-1-day ahead in the grid points and at wind farm locations

WSin = load('GridData.mat');
WSoos = load('GridData_oos.mat');
datesin = WSin.dat;
valuesin = WSin.val;
datesoos = WSoos.dat;
valuesoos = WSoos.val;

lats = WSin.lats;
lons = WSin.lons;
load('NewCoords.mat','newlat','newlon')

%%  Load the kriging maps results
%(1) m
aux = load('KriggedResults_m.mat');         params.m = aux.krigged_m';

%(2) Alpha1, Alpha2
aux = load('KriggedResults_alpha1.mat');    params.alpha1 = aux.krigged_alpha1';
aux = load('KriggedResults_alpha2.mat');    params.alpha2 = aux.krigged_alpha2;

% (3) Beta1, Beta2,...Beta24
aux = load('KriggedResults_beta1.mat');     params.beta1 = aux.krigged_beta1';
aux = load('KriggedResults_beta2.mat');     params.beta2 = aux.krigged_beta2';
aux = load('KriggedResults_beta3.mat');     params.beta3 = aux.krigged_beta3';
aux = load('KriggedResults_beta4.mat');     params.beta4 = aux.krigged_beta4';
aux = load('KriggedResults_beta5.mat');     params.beta5 = aux.krigged_beta5';
aux = load('KriggedResults_beta6.mat');     params.beta6 = aux.krigged_beta6';
aux = load('KriggedResults_beta7.mat');     params.beta7 = aux.krigged_beta7';
aux = load('KriggedResults_beta8.mat');     params.beta8 = aux.krigged_beta8';
aux = load('KriggedResults_beta9.mat');     params.beta9 = aux.krigged_beta9';
aux = load('KriggedResults_beta10.mat');    params.beta10 = aux.krigged_beta10';
aux = load('KriggedResults_beta11.mat');    params.beta11 = aux.krigged_beta11';
aux = load('KriggedResults_beta12.mat');    params.beta12 = aux.krigged_beta12';
aux = load('KriggedResults_beta13.mat');    params.beta13 = aux.krigged_beta13';
aux = load('KriggedResults_beta14.mat');    params.beta14 = aux.krigged_beta14';
aux = load('KriggedResults_beta15.mat');    params.beta15 = aux.krigged_beta15';
aux = load('KriggedResults_beta16.mat');    params.beta16 = aux.krigged_beta16';
aux = load('KriggedResults_beta17.mat');    params.beta17 = aux.krigged_beta17';
aux = load('KriggedResults_beta18.mat');    params.beta18 = aux.krigged_beta18';
aux = load('KriggedResults_beta19.mat');    params.beta19 = aux.krigged_beta19';
aux = load('KriggedResults_beta20.mat');    params.beta20 = aux.krigged_beta20';
aux = load('KriggedResults_beta21.mat');    params.beta21 = aux.krigged_beta21';
aux = load('KriggedResults_beta22.mat');    params.beta22 = aux.krigged_beta22';
aux = load('KriggedResults_beta23.mat');    params.beta23 = aux.krigged_beta23';
aux = load('KriggedResults_beta24.mat');    params.beta24 = aux.krigged_beta24';
clear aux

%(4) A's
a = load('KriggedResults_a0.mat');      params.a0 = a.krigged_a0;
a = load('KriggedResults_a1.mat');      params.a1 = a.krigged_a1;
a = load('KriggedResults_a2.mat');      params.a2 = a.krigged_a2;
a = load('KriggedResults_a3.mat');      params.a3 = a.krigged_a3;
a = load('KriggedResults_a4.mat');      params.a4 = a.krigged_a4;
a = load('KriggedResults_a5.mat');      params.a5 = a.krigged_a5;
a = load('KriggedResults_a6.mat');      params.a6 = a.krigged_a6;
a = load('KriggedResults_a7.mat');      params.a7 = a.krigged_a7;
a = load('KriggedResults_a8.mat');      params.a8 = a.krigged_a8;
a = load('KriggedResults_a9.mat');      params.a9 = a.krigged_a9;
a = load('KriggedResults_a10.mat');     params.a10 = a.krigged_adiez;
a = load('KriggedResults_a11.mat');     params.a11 = a.krigged_aonce;
a = load('KriggedResults_a12.mat');     params.a12 = a.krigged_adoce;
clear a

%(5) p's
p = load('CokrigingResults_SeasonalityDense.mat');      params.p = p.krigged_seasonality;
clear p

%(6) delta1, gamma0, gamma1,from AR(1)-GJRGARCH(1,1)
omega = load('KriggedResults_omega.mat');               params.omega = omega.krigged_omega';
delta = load('KriggedResults_delta1.mat');              params.delta = delta.krigged_delta1';
gamma = load('KriggedResults_gamma0.mat');              params.gamma1 = gamma.krigged_gamma0';
gamma = load('KriggedResults_gamma1.mat');              params.gamma0 = gamma.krigged_gamma1';
lambda = load('KriggedResults_lambda.mat');             params.lambda = lambda.krigged_lambda;
skew = load('KriggedResults_skew.mat');                params.skew = skew.krigged_skew';
shape = load('KriggedResults_shape.mat');               params.shape = shape.krigged_shape';
clear delta gamma omega skew shape lambda

%(7) q's
%q = load('CokrigingResults_SigmaBSBDense.mat');         params.q = q.krigged_sigmaBSB;
q = load('CokrigingResults_SigmaBSBDense_NegativeValuesRemovedComplete.mat');   params.q = q.krigged_sigmaBSB;


%% Make wind speed predictions at points close to wind farms
% Point 1: (35.5,-118)
idx = find(lats==35.5 & lons==-118);
newidx = find(newlat==35.5 & newlon==-118);
ws1 = valuesin(:,idx);
ws1oos = valuesoos(:,idx);

ws1Fcast = calculate_WS_hAhead(35.5, -118, newlat, newlon, params, datesin, ws1, 0, 24, '');
%ws1Fcast = calculate_WS_fullPath(35.5, -118, newlat, newlon, params, datesin, ws1, 1, 24, '');
ws1Fcastoos = calculate_WS(35.5, -118, newlat, newlon, params, datesoos, ws1oos, 0, 24, '');



%Point 2: (35, -118)
idx = find(lats==35 & lons==-118);
ws2 = valuesin(:,idx);
ws2oos = valuesoos(:,idx);

ws2Fcast = calculate_WS_hAhead(35, -118, newlat, newlon, params, datesin, ws2, 0, 24, '');
ws2Fcastoos = calculate_WS_hAhead(35, -118, newlat, newlon, params, datesoos, ws2oos, 0, 24, '');

%Point 3: (35, -118.5)
idx = find(lats==35 & lons==-118.5);
newidx = find(newlat==35 & newlon==-118.5);
ws3 = valuesin(:,idx);
ws3oos = valuesoos(:,idx);
% Make some fixes for the convergence
params.alpha1(newidx) = 0.940439306236021;
params.alpha2(newidx) = 0;
params.beta1(newidx) = 0.307026449;            params.beta2(newidx) = 0.032967373;
params.beta3(newidx) = -0.025111638;           params.beta4(newidx) = -0.058689716;
params.beta5(newidx) = -0.045511227;           params.beta6(newidx) = -0.044392088;
params.beta7(newidx) = -0.026751659;           params.beta8(newidx) = -0.020280487;
params.beta9(newidx) = -0.021554709;           params.beta10(newidx) = -0.014643684;
params.beta11(newidx) = -0.009044698;          params.beta12(newidx) = -0.007469379;
params.beta13(newidx) = -0.01268132;           params.beta14(newidx) = -0.0051046;
params.beta15(newidx) = 0.001738729;           params.beta16(newidx) = -0.000891066;
params.beta17(newidx) = 0.002130784;           params.beta18(newidx) = -0.007523688;
params.beta19(newidx) = 0.011430292;           params.beta20(newidx) = 0.003704982;
params.beta21(newidx) = 0.010673029;           params.beta22(newidx) = 0.019114414;
params.beta23(newidx) = 0.032658824;           params.beta24(newidx) = 0.055753333;

ws3Fcast = calculate_WS_hAhead(35, -118.5, newlat, newlon, params, datesin, ws3, 0, 24, '');
ws3Fcastoos = calculate_WS_hAhead(35, -118.5, newlat, newlon, params, datesoos, ws3oos, 0, 24, '');


%Point 4: (34,-116.5)
idx = find(lats==34 & lons==-116.5);
ws4 = valuesin(:,idx);
ws4oos = valuesoos(:,idx);

ws4Fcast = calculate_WS_hAhead(34, -116.5, newlat, newlon, params, datesin, ws4, 0, 24, '');
ws4Fcastoos = calculate_WS_hAhead(34, -116.5, newlat, newlon, params, datesoos, ws4oos, 0, 24, '');

%Point 5: (33,-115.5)
idx = find(lats==33 & lons==-115.5);
ws5 = valuesin(:,idx);
ws5oos = valuesoos(:,idx);

ws5Fcast = calculate_WS_hAhead(33, -115.5, newlat, newlon, params, datesin, ws5, 0, 24, '');
ws5Fcastoos = calculate_WS_hAhead(33, -115.5, newlat, newlon, params, datesoos, ws5oos, 0, 24, '');

%Point 6: (33, -116.5)
idx = find(lats==33 & lons==-116.5);
ws6 = valuesin(:,idx);
ws6oos = valuesoos(:,idx);

ws6Fcast = calculate_WS_hAhead(33, -116.5, newlat, newlon, params, datesin, ws6, 0, 24, '');
ws6Fcastoos = calculate_WS_hAhead(33, -116.5, newlat, newlon, params, datesoos, ws6oos, 0, 24, '');

%% Make wind speed predictions at wind farm locations

% Wind Farm 1: Sky River Project + NorthSkyRiver=370 MW
lat1 = 35.2864;
lon1 = -118.1948;
cap1 = [370 370 370 370 370 370];

%
idxa = find(lats==35 & lons==-118);
wsa = valuesin(:,idxa);
idxb = find(lats==35.5 & lons==-118);
wsb = valuesin(:,idxb);
idxc = find(lats==35.5 & lons==-118.5);
wsc = valuesin(:,idxc);
idxd = find(lats==35 & lons==-118.5);
wsd = valuesin(:,idxb);

ws = 0.25.*wsa+0.25.*wsb+0.25.*wsc+0.25.*wsd;
%
%[wf1, ~,~,m] = calculate_WS(lat1, lon1, newlat, newlon, params, datesin, ws1, 0, 6, '');
 wf1 = calculate_WS_hAhead(lat1, lon1, newlat, newlon, params, datesin, ws1, 0, 24, '');
 wf1Test = calculate_WS_hAhead(lat1, lon1, newlat, newlon, params, datesin, ws, 0, 24, '');
 [newW,Simepsilon,Varepsilon,m, ~] = calculate_WSPath3(35.2864,  -118.1948, newlat, newlon, params, datesin, ws1, 1,1,1,'');
%[wf1oos, ~,~,m] = calculate_WS(lat1, lon1, newlat, newlon, params, datesoos, ws1oos, 0, 6, '');
wf1oos = calculate_WS_hAhead(lat1, lon1, newlat, newlon, params, datesoos, ws1oos, 0, 24, '');

% Wind Farm 2: Alta + Oasis = 2490 MW
lat2 = [35.0634 35.0619 35.0069];
lon2 = [-118.3742 -118.2931 -118.2246];
cap2 = [2094 2094 2094 2287 2418 2490];

% [wf21, ~,~,m] = calculate_WS(lat2(1), lon2(1), newlat, newlon, params, datesin, ws2, 0, 6, '');
wf21 = calculate_WS_hAhead(lat2(1), lon2(1), newlat, newlon, params, datesin, ws2, 0, 24, '');
% [wf21oos, ~,~,m] = calculate_WS(lat2(1), lon2(1), newlat, newlon, params, datesoos, ws2oos, 0, 6, '');
wf21oos = calculate_WS_hAhead(lat2(1), lon2(1), newlat, newlon, params, datesoos, ws2oos, 0, 24, '');

% [wf22, ~,~,m] = calculate_WS(lat2(2), lon2(2), newlat, newlon, params, datesin, ws2, 0, 6, '');
 wf22 = calculate_WS_hAhead(lat2(2), lon2(2), newlat, newlon, params, datesin, ws2, 0, 24, '');
% [wf22oos, ~,~,m] = calculate_WS(lat2(2), lon2(2), newlat, newlon, params, datesoos, ws2oos, 0, 6, '');
 wf22oos = calculate_WS_hAhead(lat2(2), lon2(2), newlat, newlon, params, datesoos, ws2oos, 0, 24, '');

% [wf23, ~,~,m] = calculate_WS(lat2(3), lon2(3), newlat, newlon, params, datesin, ws2, 0, 6, '');
 wf23 = calculate_WS_hAhead(lat2(3), lon2(3), newlat, newlon, params, datesin, ws2, 0, 24, '');
% [wf23oos, ~,~,m] = calculate_WS(lat2(3), lon2(3), newlat, newlon, params, datesoos, ws2oos, 0, 6, '');
 wf23oos = calculate_WS_hAhead(lat2(3), lon2(3), newlat, newlon, params, datesoos, ws2oos, 0, 24, '');

% Wind Farm 3: Pacific Wind + Manzana = 333 MW
lat3 = 34.9095;
lon3 = -118.4410;
cap3 = [333 333 333 333 333 333];

%[wf3, ~,~,m] = calculate_WS(lat3, lon3, newlat, newlon, params, datesin, ws3, 0, 6, '');
 wf3 = calculate_WS_hAhead(lat3, lon3, newlat, newlon, params, datesin, ws3, 0, 24, '');
%[wf3oos, ~,~,m] = calculate_WS(lat3, lon3, newlat, newlon, params, datesoos, ws3oos, 0, 6, '');
 wf3oos = calculate_WS_hAhead(lat3, lon3, newlat, newlon, params, datesoos, ws3oos, 0, 24, '');

% Wind Farm 4: Mountain View + Palm Springs Repower = 648 MW
lat4 = [33.9138];
lon4 = [-116.5871];
cap4 = [606 608 608 609 609 648];

%[wf4, ~,~,m] = calculate_WS(lat4, lon4, newlat, newlon, params, datesin, ws4, 0, 6, '');
 wf4 = calculate_WS_hAhead(lat4, lon4, newlat, newlon, params, datesin, ws4, 0, 24, '');
%[wf4oos, ~,~,m] = calculate_WS(lat4, lon4, newlat, newlon, params, datesoos, ws4oos, 0, 6, '');
 wf4oos = calculate_WS(lat4, lon4, newlat, newlon, params, datesoos, ws4oos, 0, 24, '');

% Wind Farm 5: Ocotillo = 264 MW
lat5 = [32.7478];
lon5 = [-116.0717];
cap5 = [264 264 264 264 264 264];

%[wf5, ~,~,m] = calculate_WS(lat5, lon5, newlat, newlon, params, datesin, ws5, 0, 6, '');
 wf5 = calculate_WS_hAhead(lat5, lon5, newlat, newlon, params, datesin, ws5, 0, 24, '');
%[wf5oos, ~,~,m] = calculate_WS(lat5, lon5, newlat, newlon, params, datesoos, ws5oos, 0, 6, '');
 wf5oos = calculate_WS_hAhead(lat5, lon5, newlat, newlon, params, datesoos, ws5oos, 0, 24, '');

% Wind Farm 6: Tule = 183 MW
lat6 = [32.7471];
lon6 = [-116.2829];
cap6 = [52 52 52 183 183 183];

%[wf6, ~,~,m] = calculate_WS(lat6, lon6, newlat, newlon, params, datesin, ws6, 0, 6, '');
 wf6 = calculate_WS_hAhead(lat6, lon6, newlat, newlon, params, datesin, ws6, 0, 24, '');
%[wf6oos, ~,~,m] = calculate_WS(lat6, lon6, newlat, newlon, params, datesoos, ws6oos, 0, 6, '');
 wf6oos = calculate_WS_hAhead(lat6, lon6, newlat, newlon, params, datesoos, ws6oos, 0, 24, '');


%% Save the forecasted wind speeds at wind farm locations
save('Actualh__GridPoint_WSCAISO_insample.mat','ws1','ws2','ws3','ws4','ws5','ws6','cap1','cap2','cap3','cap4','cap5','cap6','datesin')
save('Actualh_GridPoint_WSCAISO_oos.mat','ws1oos','ws2oos','ws3oos','ws4oos','ws5oos','ws6oos','cap1','cap2','cap3','cap4','cap5','cap6','datesoos')

save('Forecasted_24h_GridPoint_WSCAISO_insample.mat','ws1Fcast','ws2Fcast','ws3Fcast','ws4Fcast','ws5Fcast','ws6Fcast','cap1','cap2','cap3','cap4','cap5','cap6','datesin')
save('Forecasted_24h_GridPoint_WSCAISO_oos.mat','ws1Fcastoos','ws2Fcastoos','ws3Fcastoos','ws4Fcastoos','ws5Fcastoos','ws6Fcastoos','cap1','cap2','cap3','cap4','cap5','cap6','datesoos')

save('Forecasted_24h_WSCAISO_insample.mat','wf1','wf21','wf22','wf23','wf3','wf4','wf5','wf6','cap1','cap2','cap3','cap4','cap5','cap6','datesin')
save('Forecasted_24h_WSCAISO_oos.mat','wf1oos','wf21oos','wf22oos','wf23oos','wf3oos','wf4oos','wf5oos','wf6oos','cap1','cap2','cap3','cap4','cap5','cap6','datesoos')


