%% Start Analyzing the hourly time series
load('GridData.mat')



v = val;
d = dat;


%Do the Weibull fit
for st=1:85
    data = double(v(:,st));
    [par,~]= mle(data,'distribution','weibull','alpha',0.05);
    c(st) = par(1); %scale param
    k(st) = par(2); %shape param
    % Plot
    x = min(data):range(data)/100:max(data);
    [bincount,binpos] = hist(data,min(100,numel(data)/5));
    y = pdf('weibull',x,c(st),k(st));
    figure; hold on
    % binwidth = binpos(2) - binpos(1);
    % histarea = binwidth*sum(bincount);
    bincount= bincount/trapz(binpos,bincount); % scaled frequencies
    data= bar(binpos,bincount,'FaceColor',[.8 .8 .8],'EdgeColor',[0.71 0.71 0.71],'BarWidth',1);
    model= plot(x,y,'k','LineWidth',2);
    xlabel('Data'); ylabel('PDF')
    legend([data,model],'Data','Weibull'); legend('boxoff');
    title(sprintf('Distribution in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
end

%% Transform wind speeds by the Weibull shape param/3.6 (Nfaoui, Torres etc.)
m = k'./3.6;
for st=1:85
    newv(:,st) = v(:,st).^m(st);
end

%Do the Normal fit
for st=1:85
    data = double(newv(:,st));
    [par,~]= mle(data,'distribution','normal','alpha',0.05);
    c(st) = par(1); %mean param
    k(st) = par(2); %var param
    % Plot
    x = min(data):range(data)/100:max(data);
    [bincount,binpos] = hist(data,min(100,numel(data)/5));
    y = pdf('normal',x,c(st),k(st));
    figure; hold on
    % binwidth = binpos(2) - binpos(1);
    % histarea = binwidth*sum(bincount);
    bincount= bincount/trapz(binpos,bincount); % scaled frequencies
    data= bar(binpos,bincount,'FaceColor',[.8 .8 .8],'EdgeColor',[0.71 0.71 0.71],'BarWidth',1);
    model= plot(x,y,'k','LineWidth',2);
    xlabel('Data'); ylabel('PDF')
    legend([data,model],'Data','Normal'); legend('boxoff');
    title(sprintf('New distribution in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
end


% val = v; save('Hourly_transformed.mat','dat','val')
%% Deseasonalized
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

 
numCycles2 = 2;
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

% A = fitlm(eff(:,2:end),v);

% R = regstats(v(:,st),eff(:,2:end));

% sprintf('AIC for %d yearly and %d daily is: %.4f',numCycles, numCycles2,aic(numCycles, numCycles2))
% n = length(v);
% SSE = R.fstat.sse;
% RSS = sum(R.r.^2);
% p = size(eff,2);
% AIC(st) = log(SSE/n) + (n+2*p)/n;
% AIC(st) = n*log(RSS/n)+2*p;
% BIC(st, numCycles) = log(SSE/n) + p*log(n)/n;
% Rsq(st, numCycles) = R.rsquare;
% MSE(st, numCycles) = R.mse;
% end


% m = [];
% for st=1:85
%     m(st) = find(AIC(st,:)==min(AIC(st,:)));
% end

% load('Daily_stations_seasonality.mat')
%%
% [b, fitinfo] = lasso(eff, v);
sest = b(1).*ones(size(vv));
for i=1:size(eff,2)-1
    sest = sest+b(i+1).*eff(:,i+1);
end

figure;tSeries(d,v(:,st)).plot('Color',[0.005 0.005 0.005]);
hold on;
plot(d,sest,'r', 'LineWidth',3);
xlabel('Date')
ylabel('Hourly seasonal effect')
title(sprintf('Seasonality in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))


allseasonalParams = [allseasonalParams b];
% alpha.(char(field)) = b;
% Gamma.(char(field)) = sest;
% Stats.(char(field)) = R;
ws(:,st) = v(:,st)-sest;
clear sest R cycles cyclesH eff b
end
val = ws;
save('Hourly_deseasonalized.mat','dat','val','lats', 'lons')
%% Plot distribution of the deseasonalized data
for st=1:85
    data = double(val(:,st));
    [par,~]= mle(data,'distribution','normal','alpha',0.05);
    mu(st) = par(1); %scale param
    sigma(st) = par(2); %shape param
    % Plot
    x = min(data):range(data)/100:max(data);
    [bincount,binpos] = hist(data,min(100,numel(data)/5));
    y = pdf('normal',x,mu(st),sigma(st));
    figure; hold on
    % binwidth = binpos(2) - binpos(1);
    % histarea = binwidth*sum(bincount);
    bincount= bincount/trapz(binpos,bincount); % scaled frequencies
    data= bar(binpos,bincount,'FaceColor',[.8 .8 .8],'EdgeColor',[0.71 0.71 0.71],'BarWidth',1);
    model= plot(x,y,'k','LineWidth',2);
    xlabel('Data'); ylabel('PDF')
    legend([data,model],'Data','Normal'); legend('boxoff');
    title(sprintf('Deseasonalized data in point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)))
end

%% Collapse all values in 1 day to observe daily profile
localvalid = Timezone(dat,'UTC', 'PST');
for h=0:23
    idxS = hour(localvalid)==h & ismember(month(localvalid), 3:9);
    idxW = hour(localvalid)==h & ~ismember(month(localvalid), 3:9);
    for st=1:85
        vhS(h+1,st) = mean(v(idxS,st));
        vhW(h+1,st) = mean(v(idxW,st));
    end
end


for st=1:85
figure; 
subplot(2,1,1)
plot(0:23, vhW(:,st),'*-k')
xlabel('Hour (Winter)')
ylabel('Wind speed')
subplot(2,1,2)
plot(0:23, vhS(:,st),'*-k')
xlabel('Hour (Summer)')
ylabel('Wind speed')
a = axes;
t = title(sprintf('Hourly profile in grid point %d (lat:%.2f, lon:%.2f)',st, lats(st), lons(st)));
a.Visible = 'off'; % set(a,'Visible','off');
t.Visible = 'on'; % set(t1,'Visible','on');

end