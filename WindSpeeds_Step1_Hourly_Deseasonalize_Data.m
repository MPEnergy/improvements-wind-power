%% Start Analyzing the hourly time series
clear all
load('GridData.mat')


% 1. Check ACF, PACF
v = val;
d = dat;

% Load the m exponents for transform

T = xlsread('Transformations.xlsx', 'Exponent');
m = T(:,7);

%1.transform
for st=1:85
    v(:, st) = val(:,st).^m(st);
end



%% Remove first the yearly seasonality
ws = nan(size(v));
valid = d;
vv = round(valid.*24)./24-datenum(year(valid),1,1);
vh = hour(valid);
allseasonalParams = [];
for st=1:85
numCycles = 6;
vv = round(valid.*24)./24-datenum(year(valid),1,1);
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
vv = hour(valid);
cyclesH = repmat(struct('dat',[],'name',[]),1,2*numCycles2);
for ii = 1:numCycles2
    tc = cos(2*pi*vv/24*ii);
    ts = sin(2*pi*vv/24*ii);
    cyclesH(2*(ii-1)+1).dat = tc;
    cyclesH(2*(ii-1)+1).name = ['cosH' num2str(ii)];
    cyclesH(2*ii).dat = ts;
    cyclesH(2*ii).name = ['sinH' num2str(ii)];
end

eff = [ones(size(valid)) [cycles.dat cyclesH.dat]];

field = sprintf('st%d',st);
b = regress(v(:,st), eff);
%%

sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end

ws(:,st) = v(:,st)-sest;

allseasonalParams = [allseasonalParams b];

figure;tSeries(d,v(:,st)).plot('Color',[0.005 0.005 0.005]);
hold on;
plot(d,sest,'r', 'LineWidth',3);
xlabel('Date')
ylabel('Hourly seasonal effect')
title(sprintf('Seasonality in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))


end

%% Build monthly-specific hourly profiles

%localvalid = Timezone(valid,'UTC', 'PST');
localvalid= valid;
v = ws;
for h=0:23
    idxJan = hour(localvalid)==h & month(localvalid)==1;
    idxFeb = hour(localvalid)==h & month(localvalid)==2;
    idxMar = hour(localvalid)==h & month(localvalid)==3;
    idxApr = hour(localvalid)==h & month(localvalid)==4;
    idxMay = hour(localvalid)==h & month(localvalid)==5;
    idxJun = hour(localvalid)==h & month(localvalid)==6;
    idxJul = hour(localvalid)==h & month(localvalid)==7;
    idxAug = hour(localvalid)==h & month(localvalid)==8;
    idxSep = hour(localvalid)==h & month(localvalid)==9;
    idxOct = hour(localvalid)==h & month(localvalid)==10;
    idxNov = hour(localvalid)==h & month(localvalid)==11;
    idxDec = hour(localvalid)==h & month(localvalid)==12;
    
    for st=1:85
        vhJan(h+1,st) = mean(v(idxJan,st));
        vhFeb(h+1,st) = mean(v(idxFeb,st));
        vhMar(h+1,st) = mean(v(idxMar,st));
        vhApr(h+1,st) = mean(v(idxApr,st));
        vhMay(h+1,st) = mean(v(idxMay,st));
        vhJun(h+1,st) = mean(v(idxJun,st));
        vhJul(h+1,st) = mean(v(idxJul,st));
        vhAug(h+1,st) = mean(v(idxAug,st));
        vhSep(h+1,st) = mean(v(idxSep,st));
        vhOct(h+1,st) = mean(v(idxOct,st));
        vhNov(h+1,st) = mean(v(idxNov,st));
        vhDec(h+1,st) = mean(v(idxDec,st));
    end
end

% Make some charts
for st=1:85
    figure
    subplot(3,2,1)
    plot(1:24, vhJan(:,st),'-*k')
    title('January')
    subplot(3,2,2)
    plot(1:24, vhFeb(:,st),'-*k')
    title('February')
    subplot(3,2,3)
    plot(1:24, vhMar(:,st),'-*k')
    title('March')
    subplot(3,2,4)
    plot(1:24, vhApr(:,st),'-*k')
    title('April')
    subplot(3,2,5)
    plot(1:24, vhMay(:,st),'-*k')
    title('May')
    subplot(3,2,6)
    plot(1:24, vhJun(:,st),'-*k')
    title('June')
    
    figure;
    subplot(3,2,1)
    plot(1:24, vhJul(:,st),'-*k')
    title('July')
    subplot(3,2,2)
    plot(1:24, vhAug(:,st),'-*k')
    title('August')
    subplot(3,2,3)
    plot(1:24, vhSep(:,st),'-*k')
    title('September')
    subplot(3,2,4)
    plot(1:24, vhOct(:,st),'-*k')
    title('October')
    subplot(3,2,5)
    plot(1:24, vhNov(:,st),'-*k')
    title('November')
    subplot(3,2,6)
    plot(1:24, vhDec(:,st),'-*k')
    title('December')
end
% Build synthetic daily time series
syn = nan(length(localvalid), 85);
for st=1:85
    for h=0:23
        idxJan = hour(localvalid)==h & month(localvalid)==1;
        idxFeb = hour(localvalid)==h & month(localvalid)==2;
        idxMar = hour(localvalid)==h & month(localvalid)==3;
        idxApr = hour(localvalid)==h & month(localvalid)==4;
        idxMay = hour(localvalid)==h & month(localvalid)==5;
        idxJun = hour(localvalid)==h & month(localvalid)==6;
        idxJul = hour(localvalid)==h & month(localvalid)==7;
        idxAug = hour(localvalid)==h & month(localvalid)==8;
        idxSep = hour(localvalid)==h & month(localvalid)==9;
        idxOct = hour(localvalid)==h & month(localvalid)==10;
        idxNov = hour(localvalid)==h & month(localvalid)==11;
        idxDec = hour(localvalid)==h & month(localvalid)==12;
        syn(idxJan,st) = repmat(vhJan(h+1,st), length(find(idxJan)), 1);
        syn(idxFeb,st) = repmat(vhFeb(h+1,st), length(find(idxFeb)), 1);
        syn(idxMar,st) = repmat(vhMar(h+1,st), length(find(idxMar)), 1);
        syn(idxApr,st) = repmat(vhApr(h+1,st), length(find(idxApr)), 1);
        syn(idxMay,st) = repmat(vhMay(h+1,st), length(find(idxMay)), 1);
        syn(idxJun,st) = repmat(vhJun(h+1,st), length(find(idxJun)), 1);
        syn(idxJul,st) = repmat(vhJul(h+1,st), length(find(idxJul)), 1);
        syn(idxAug,st) = repmat(vhAug(h+1,st), length(find(idxAug)), 1);
        syn(idxSep,st) = repmat(vhSep(h+1,st), length(find(idxSep)), 1);
        syn(idxOct,st) = repmat(vhOct(h+1,st), length(find(idxOct)), 1);
        syn(idxNov,st) = repmat(vhNov(h+1,st), length(find(idxNov)), 1);
        syn(idxDec,st) = repmat(vhDec(h+1,st), length(find(idxDec)), 1);
    end
end




st = 1;
figure
plot(1:24, vhJan(:, st), '-*')
title(sprintf('January profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
figure
plot(1:24, vhFeb(:, st), '-*')
title(sprintf('February profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
figure
plot(1:24, vhMar(:, st), '-*')
title(sprintf('March profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
figure
plot(1:24, vhApr(:, st), '-*')
title(sprintf('April profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhMay(:, st), '-*')
title(sprintf('May profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhJun(:, st), '-*')
title(sprintf('June profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhJul(:, st), '-*')
title(sprintf('July profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhAug(:, st), '-*')
title(sprintf('August profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhSep(:, st), '-*')
title(sprintf('September profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhOct(:, st), '-*')
title(sprintf('October profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhNov(:, st), '-*')
title(sprintf('November profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

figure
plot(1:24, vhDec(:, st), '-*')
title(sprintf('December profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

%% Make some plots
st = 66;
figure; hold on;
tSeries(localvalid, v(:,st)).plot('k');
tSeries(localvalid, syn(:,st)).plot('r');
title(sprintf('Seasonality in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))


%% Final deseasonalization
val = v-syn;
save('Hourly_deseasonalized_test.mat','dat','val')
