%% Model the residuals after ARMA: 
% epsilon=sigma*xi
% sigma^2=sigma_BSB^2*sigma_GARCH^2
% sigma_BSB^2 Fourier
% sigma_GARCH^2 GARCH(1,1)
clear all
load('Hourly_deseasonalized.mat', 'dat', 'lats', 'lons');
p = xlsread('Residuals_ARMA(2,24).xlsx', 'CorrectResiduals');

Res = p(1:end, 2:end);


%val=Res;save('Hourly_residuals_after_ARMA.mat','dat','val')

%% Seasonal part


m = month(dat);

for st =1:85
    monthlyAverage(:, st) = accumarray(m, Res(:,st).^2, [], @mean);
    figure;
    plot(1:12, monthlyAverage(:,st), '-*k')
end

% Fit with Fourier series
numCycles = 1;
vv = [1:12]';
for ii = 1:numCycles
    tc = cos(2*pi*vv/12*ii);
    ts = sin(2*pi*vv/12*ii);
    cycles(:, ii) = tc';
    cycles(:, ii+1) = ts';

end

eff = [ones(12,1) cycles];
regParams = [];
for st=1:85
b = regress(monthlyAverage(:,st), eff);
regParams = [regParams, b];
sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end


figure;
hold on;
plot(1:12, monthlyAverage(:,st), '-*k')
plot(1:12, sest, 'r')
xlabel('Week no')
ylabel('Variance') 
title(sprintf('Fitted variance in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))





sigma2(:, st) = sest;

end

%% Extend this over the whole in-sample period
for st=1:85
    epsilon(:,st) = Res(:,st)./sqrt(sigma2(m,st));
end

%%  Weekly
week = WeekNumber(dat);

for st =1:85
    weeklyAverage(:, st) = accumarray(week, Res(:,st).^2, [], @mean);
    figure;
    plot(1:53, weeklyAverage(:,st), '-*k')
end

% Fit with Fourier series
numCycles = 2;
vv = [1:53]';
for ii = 1:numCycles
    tc = cos(2*pi*vv/52.18*ii);
    ts = sin(2*pi*vv/52.18*ii);
    cycles(:, ii) = tc';
    cycles(:, ii+1) = ts';

end

eff = [ones(53,1) cycles];
regParams = [];
for st=1:85
b = regress(weeklyAverage(:,st), eff);
regParams = [regParams, b];
sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end

figure;
hold on;
plot(1:53, weeklyAverage(:,st), '-*k')
plot(1:53, sest, 'r')
xlabel('Week no')
ylabel('Variance') 
title(sprintf('Fitted variance in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

sigma2(:, st) = sest;

end

%% Extend this over the whole in-sample period
for st=1:85
    epsilon(:,st) = Res(:,st)./sqrt(sigma2(week,st));
end

val = epsilon; save('Hourly_residuals_weeklyseasonality.mat','dat','val')


%%  Daily
doy = DayOfYear(floor(dat));

for st =1:85
    dailyAverage(:, st) = accumarray(doy, Res(:,st).^2, [], @mean);
    figure;
    plot(1:366, dailyAverage(:,st), '-*k')
end

% Fit with Fourier series
numCycles = 2;
vv = [1:366]';
for ii = 1:numCycles
    tc = cos(2*pi*vv/365.25*ii);
    ts = sin(2*pi*vv/365.25*ii);
    cycles(:, ii) = tc';
    cycles(:, ii+1) = ts';

end

eff = [ones(366,1) cycles];
regParams = [];
for st=1:85
b = regress(dailyAverage(:,st), eff);
regParams = [regParams, b];
sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end

figure;
hold on;
plot(1:366, dailyAverage(:,st), '-*k')
plot(1:366, sest, 'r')
xlabel('Day no')
ylabel('Variance') 
title(sprintf('Fitted variance in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

sigma2(:, st) = sest;

end

%% Extend this over the whole in-sample period
for st=1:85
    epsilon(:,st) = Res(:,st)./sqrt(sigma2(doy,st));
end

%% Hourly
newdat = Timezone(dat, 'UTC', 'PST');
hod = hour(dat)+ones(size(dat));

for st =1:85
    hourlyAverage(:, st) = accumarray(hod, Res(:,st).^2, [], @mean);
    figure;
    plot(0:23, hourlyAverage(:,st), '-*k')
    xlabel('Hour')
    ylabel('Variance')
    %title(sprintf(' in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
end

% Fit with Fourier series
numCycles = 2;
vv = [1:24]';
for ii = 1:numCycles
    tc = cos(2*pi*vv/24*ii);
    ts = sin(2*pi*vv/24*ii);
    cycles(:, ii) = tc';
    cycles(:, ii+1) = ts';

end

eff = [ones(24,1) cycles];
regParams = [];
for st=1:85
b = regress(hourlyAverage(:,st), eff);
regParams = [regParams, b];
sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end

figure;
hold on;
plot(0:23, hourlyAverage(:,st), '-*k')
plot(0:23, sest, 'r')
xlabel('Hour')
ylabel('Variance') 
title(sprintf('Fitted variance in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))

sigma2(:, st) = sest;

end


%% Extend this over the whole in-sample period
for st=1:85
    epsilon(:,st) = Res(:,st)./sqrt(sigma2(hod,st));
end




val = epsilon; save('Hourly_residuals_hourlyseas.mat','dat','val')

%% Check hourly bias per month
v = Res.^2;
localvalid= dat;

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
    
    %sgtitle(sprintf('Point %d', st))
    
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
    %sgtitle(sprintf('Point %d', st))
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

for st=1:85
    s(:,st) = Res(:,st)./sqrt(syn(:,st));
end

val = s;
save('Hourly_residuals_sigmaGARCH.mat','dat','val')

%% Each hour in 1 year

for st=1:85
    vhJan(:,st) = hpfilter(vhJan(:,st),1);
    vhFeb(:,st) = hpfilter(vhFeb(:,st),1);
    vhMar(:,st) = hpfilter(vhMar(:,st),1);
    vhApr(:,st) = hpfilter(vhApr(:,st),1);
    vhMay(:,st) = hpfilter(vhMay(:,st),1);
    vhJun(:,st) = hpfilter(vhJun(:,st),1);
    vhJul(:,st) = hpfilter(vhJul(:,st),1);
    vhAug(:,st) = hpfilter(vhAug(:,st),1);
    vhSep(:,st) = hpfilter(vhSep(:,st),1);
    vhOct(:,st) = hpfilter(vhOct(:,st),1);
    vhNov(:,st) = hpfilter(vhNov(:,st),1);
    vhDec(:,st) = hpfilter(vhDec(:,st),1);
end
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

for st=1:85
    s2(:,st) = Res(:,st)./sqrt(syn(:,st));
end
%% Make some plots
st = 1;
figure; hold on;
tSeries(localvalid, v(:,st)).plot('k');
tSeries(localvalid, syn(:,st)).plot('r');
title(sprintf('Seasonality in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))