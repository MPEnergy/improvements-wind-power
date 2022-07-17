%% Step2: Make spatial predictions at wind farm locations
clear all
WSin = load('GridData.mat');
WSoos = load('GridData_oos.mat');
datesin = WSin.dat;
valuesin = WSin.val;
datesoos = WSoos.dat;
valuesoos = WSoos.val;

lats = WSin.lats;
lons = WSin.lons;
load('NewCoords.mat','newlat','newlon')

% load('Hourly_deseasonalized_test.mat')
% datesin = dat;
% valuesin = val;

%% Divide the data in summer(May-September) day/night and winter(November-April) day/night
for h=0:23
    
%     fsummer = char(sprintf('summerH%d',h));
%     fwinter = char(sprintf('winterH%d',h));
%     
%     idxsummer = find(ismember(month(datesin),[5:9]) & hour(datesin)==h);
%     idxwinter = find(ismember(month(datesin),union([1:4],[10:12])) & hour(datesin)==h);
%     
%     valsin.(fsummer) = valuesin(idxsummer,:);       datin.(fsummer) = datesin(idxsummer);
%     valsin.(fwinter) = valuesin(idxwinter,:);       datin.(fwinter) = datesin(idxwinter);
    
    fJan = char(sprintf('JanH%d',h));
    fFeb = char(sprintf('FebH%d',h));
    fMar = char(sprintf('MarH%d',h));
    fApr = char(sprintf('AprH%d',h));
    fMay = char(sprintf('MayH%d',h));
    fJun = char(sprintf('JunH%d',h));
    fJul = char(sprintf('JulH%d',h));
    fAug = char(sprintf('AugH%d',h));
    fSep = char(sprintf('SepH%d',h));
    fOct = char(sprintf('OctH%d',h));
    fNov = char(sprintf('NovH%d',h));
    fDec = char(sprintf('DecH%d',h));
%     %In-sample
    idxJan = find(month(datesin)==1 & hour(datesin)==h);
    idxFeb = find(month(datesin)==2 & hour(datesin)==h);
    idxMar = find(month(datesin)==3 & hour(datesin)==h);
    idxApr = find(month(datesin)==4 & hour(datesin)==h);
    idxMay = find(month(datesin)==5 & hour(datesin)==h);
    idxJun = find(month(datesin)==6 & hour(datesin)==h);
    idxJul = find(month(datesin)==7 & hour(datesin)==h);
    idxAug = find(month(datesin)==8 & hour(datesin)==h);
    idxSep = find(month(datesin)==9 & hour(datesin)==h);
    idxOct = find(month(datesin)==10 & hour(datesin)==h);
    idxNov = find(month(datesin)==11 & hour(datesin)==h);
    idxDec = find(month(datesin)==12 & hour(datesin)==h);
    
    valsin.(fJan) = valuesin(idxJan,:);
    valsin.(fFeb) = valuesin(idxFeb,:);
    valsin.(fMar) = valuesin(idxMar,:);
    valsin.(fApr) = valuesin(idxApr,:);
    valsin.(fMay) = valuesin(idxMay,:);
    valsin.(fJun) = valuesin(idxJun,:);
    valsin.(fJul) = valuesin(idxJul,:);
    valsin.(fAug) = valuesin(idxAug,:);
    valsin.(fSep) = valuesin(idxSep,:);
    valsin.(fOct) = valuesin(idxOct,:);
    valsin.(fNov) = valuesin(idxNov,:);
    valsin.(fDec) = valuesin(idxDec,:);

    % Out-of-sample
    
%     idxsummer = find(ismember(month(datesoos),[5:9]) & hour(datesoos)==h);
%     idxwinter = find(ismember(month(datesoos),union([1:4],[10:12])) & hour(datesoos)==h);
%     
%     valsoos.(fsummer) = valuesoos(idxsummer,:);
%     valsoos.(fwinter) = valuesoos(idxwinter,:);
    
    
    
    idxJan = find(month(datesoos)==1 & hour(datesoos)==h);
    idxFeb = find(month(datesoos)==2 & hour(datesoos)==h);
    idxMar = find(month(datesoos)==3 & hour(datesoos)==h);
    idxApr = find(month(datesoos)==4 & hour(datesoos)==h);
    idxMay = find(month(datesoos)==5 & hour(datesoos)==h);
    idxJun = find(month(datesoos)==6 & hour(datesoos)==h);
    idxJul = find(month(datesoos)==7 & hour(datesoos)==h);
    idxAug = find(month(datesoos)==8 & hour(datesoos)==h);
    idxSep = find(month(datesoos)==9 & hour(datesoos)==h);
    idxOct = find(month(datesoos)==10 & hour(datesoos)==h);
    idxNov = find(month(datesoos)==11 & hour(datesoos)==h);
    idxDec = find(month(datesoos)==12 & hour(datesoos)==h);
    
    valsoos.(fJan) = valuesoos(idxJan,:);
    valsoos.(fFeb) = valuesoos(idxFeb,:);
    valsoos.(fMar) = valuesoos(idxMar,:);
    valsoos.(fApr) = valuesoos(idxApr,:);
    valsoos.(fMay) = valuesoos(idxMay,:);
    valsoos.(fJun) = valuesoos(idxJun,:);
    valsoos.(fJul) = valuesoos(idxJul,:);
    valsoos.(fAug) = valuesoos(idxAug,:);
    valsoos.(fSep) = valuesoos(idxSep,:);
    valsoos.(fOct) = valuesoos(idxOct,:);
    valsoos.(fNov) = valuesoos(idxNov,:);
    valsoos.(fDec) = valuesoos(idxDec,:);
end


%% Take 1 group of data (in-sample)/Repeat this step for out-of-sample (line116)
% % A = [];
% % C = [];
% % nugget = [];
% % for h=0:23
% %     f = char(sprintf('NovH%d',h));
% %     A = [];
% %     C = [];
% % for v=1:size(valsin.(f),1)
% % ws = valsin.(f)(v,:);
% % % s = regstats(ws',[lons lats lons.*lats lons.^2 lats.^2], 'linear');
% % % ws = s.r;
% % %% Make the semivariogram
% % for i= 1:85
% %     for j = 1:85
% %         [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
% %         B(i,j) = bearing([lats(i) lons(i)], [lats(j) lons(j)]);
% %     end 
% % end
% % 
% % Ns = unique(D);
% % edges = 0:11:800;
% % [Nb, BIN] = histc(Ns, edges);
% % 
% % % Calculate the semivariogram
% % %[g0, cf0] = calculate_semivariogram(ws, edges, D, B, 70);
% % [g0, cf0] = calculate_semivariogram(ws, edges, D);
% % % % 
% % % % 
% % % % 
% % % % %% Fit the parameters
% % % % 
% % % % % load('Semivariograms_b_coeffs.mat','edges','g0','g1','g2','g3','g4','cf0','cf1','cf2','cf3','cf4')
% % % % % A0
% % xdata = edges';
% % ydata = g0';
% % fun1 = @(x)ssevalexp(x,xdata,ydata);
% % % % 
% % % % 
% % % % % x0 = [400, 0.0015]; %% [400,10, 50, -20]
% % x0 = [400,8];
% % options = optimset('MaxFunEvals',1000000);
% % bestx1 = fminsearch(fun1,x0, options);
% % % sph = @(x,a,c0,c) c0*ones(size(x))+c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+(c0+c).*(x>a)+c0.*(x==0);
% % sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
% % gaus = @(x,a,c) c.*(1-exp(-x.^2/a.^2));
% % ex = @(x,a,c) c.*(1-exp(-x./a));
% % cub = @(x,a,c) c.*(7.*(x./a).^2-(35/4).*(x./a).^3+(7/2).*(x./a).^5-(3/4).*(x./a).^7).*(x<=a)+c.*(x>a);
% % sinehole = @(x,a,c) c.*(1-sin(pi.*x./a)./(pi.*x./a)).*(x~=0);
% % power = @(x,a,c) c.*x.^a;  %a<2, a>=0
% % % % 
% % % % %Check the fit
% % a1 = bestx1(1);
% % c1 = bestx1(2);
% % % a2 = bestx1(3);
% % % c2 = bestx1(4);
% % % % 
% % % % 
% % % yfit = sph(xdata, a1, c1);
% % % figure;
% % % plot(xdata,ydata,'*');
% % % hold on
% % % plot(xdata,yfit,'r');
% % % xlabel('distance (km)')
% % % ylabel('Semivariogram (m/s)')
% % % title(sprintf('Sine Hole Effect and Exponential function'))
% % % legend('b_0 semivariogram','Fitted semivariogram')
% % % hold off
% % % % 
% % % % 
% % % fitsemi = @(h,a1,c1) (sph(h,a1,c1));  %the fitted semivariogram function
% % % yfit = fitsemi(xdata,a1,c1);
% % % figure;
% % % 
% % % subplot(2,1,1)
% % % plot(xdata,ydata,'*');
% % % hold on
% % % plot(xdata,yfit,'r');
% % % xlabel('distance (km)')
% % % ylabel('Semivariogram')
% % % title('Fit of spatial variation for $\beta_{22}$', 'Interpreter', 'latex')
% % % legend('\beta_{22} semivariogram function','Fitted semivariogram')
% % % hold off
% % % 
% % % 
% % % xdata = edges';
% % % ydata = cf0';
% % % fitcov = @(h,a1,c1) (cf0(1) -fitsemi(h,a1,c1));
% % % 
% % % yfit = fitcov(xdata,a1,c1);
% % % subplot(2,1,2)
% % % plot(xdata,ydata,'*');
% % % hold on
% % % plot(xdata,yfit,'r');
% % % xlabel('distance (km)')
% % % ylabel('Covariance')
% % % legend('\beta_{22} covariance function','Fitted covariance')
% % % hold off
% % % % % nuggetbeta22 = 3.476;
% % % % % 
% % % % % 
% % % % 
% % % % 
% % % % 
% % % % 
% % % % %sqrt(nanmean((ydata-yfit).^2))
% % % % 
% % A = [A a1];
% % C = [C c1];
% % nugget = [nugget cf0(1)];
% % end
% % Af(h+1) = mean(A(A<800));
% % Cf(h+1) = mean(C(A<800));
% % nuggetf(h+1) = mean(nugget(A<800));
% % end

%% Make some plots
% % KrigingParams = xlsread('C:\Users\u6032456\Documents\Thesis\PhDMeeting\article3\AnisotropicSemivariograms.xlsx', 1, 'A52:D339');
% % A = KrigingParams(:,1);
% % C = KrigingParams(:,2);
% % sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);
% % 
% % figure
% % h = 0;
% % f = char(sprintf('SepH%d',h));
% % ws = valsin.(f)(99,:);
% % for i= 1:85
% %     for j = 1:85
% %         [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
% %         B(i,j) = bearing([lats(i) lons(i)], [lats(j) lons(j)]);
% %     end 
% % end
% % 
% % Ns = unique(D);
% % edges = 0:11:800;
% % [Nb, BIN] = histc(Ns, edges);
% % [g0, cf0] = calculate_semivariogram(ws, edges, D);
% % xdata = edges';
% % ydata1 = g0';
% % 
% % a1 = A(1);
% % c1 = C(1);
% % yfit = sph(xdata, a1, c1);
% % 
% % 
% % ws = valsin.(f)(120,:);
% % for i= 1:85
% %     for j = 1:85
% %         [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
% %         B(i,j) = bearing([lats(i) lons(i)], [lats(j) lons(j)]);
% %     end 
% % end
% % 
% % Ns = unique(D);
% % edges = 0:11:800;
% % [Nb, BIN] = histc(Ns, edges);
% % [g0, cf0] = calculate_semivariogram(ws, edges, D);
% % xdata = edges';
% % ydata2 = g0';
% % 
% % 
% % subplot(2,1,2)
% % hold on;
% % plot(xdata,ydata1,'*k');
% % plot(xdata,ydata2,'sk');
% % plot(xdata,yfit,'k');
% % xlabel('distance (km)')
% % ylabel('m/s')
% % title('September 15:00 PM')
% % hold off
% % 
% % 
% % h = 12;
% % f = char(sprintf('SepH%d',h));
% % ws = valsin.(f)(105,:);
% % for i= 1:85
% %     for j = 1:85
% %         [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
% %         B(i,j) = bearing([lats(i) lons(i)], [lats(j) lons(j)]);
% %     end 
% % end
% % 
% % Ns = unique(D);
% % edges = 0:11:800;
% % [Nb, BIN] = histc(Ns, edges);
% % [g0, cf0] = calculate_semivariogram(ws, edges, D);
% % xdata = edges';
% % ydata1 = g0';
% % a1 = A(13);
% % c1 = C(13);
% % yfit = sph(xdata, a1, c1);
% % 
% % 
% % ws = valsin.(f)(119,:);
% % for i= 1:85
% %     for j = 1:85
% %         [D(i,j),~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
% %         B(i,j) = bearing([lats(i) lons(i)], [lats(j) lons(j)]);
% %     end 
% % end
% % 
% % Ns = unique(D);
% % edges = 0:11:800;
% % [Nb, BIN] = histc(Ns, edges);
% % [g0, cf0] = calculate_semivariogram(ws, edges, D);
% % xdata = edges';
% % ydata2 = g0';
% % 
% % 
% % subplot(2,1,1)
% % hold on;
% % h1= plot(xdata,ydata1,'*k');
% % h2 = plot(xdata,ydata2,'sk');
% % h3 = plot(xdata,yfit,'k');
% % xlabel('distance (km)')
% % ylabel('m/s')
% % title('September 03:00 AM')
% % hold off
% % legend([h3],{'spherical semivariogram'})
% % %subtitle('Semivariogram fit for January 15:00 PM')


%% Look for trend surface


% % for h=0:23
% %     f = char(sprintf('summerH%d',h));
% %     A = [];
% %     C = [];
% % for i=1:100%length(valsin.(f))
% % ws = valsin.(f)(i,:);
% % 
% % 
% % 
% % 
% % [B,BINT,R,RINT,STATS] = regress(ws',[ones(length(ws),1) lons lats]);
% % s = regstats(ws',[lons lats lons.*lats lons.^2 lats.^2], 'linear');
% % 
% % 
% % figure(i)
% % subplot(2,1,1)
% % hist(ws)
% % subplot(2,1,2)
% % hist(s.r)
% % 
% % % figure(i)
% % % scatter(lons,lats, 60, ws, 'filled');
% % % colormap(hsv) %can be removed
% % % xlabel('Lon')
% % % ylabel('Lat')
% % % colorbar
% % 
% % 
% % 
% % end
% % end
% % 







%% Start building the kriging maps

%1. Load the parameters
months = {'Jan', 'Feb', 'Mar', 'Apr', 'May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};


KrigingParams = xlsread('AnisotropicSemivariograms.xlsx', 'PureOOS', 'A1:D289');
A = KrigingParams(:,1);
C = KrigingParams(:,2);
Nugget = KrigingParams(:,3);
sph = @(x,a,c) c.*(3.*x./(2*a)-(1/2).*(x.^3./a^3)).*(x>0 & x<=a)+c.*(x>a);


%2. Make the dense grid


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

n = length(newlat);


%3. Loop through all
% load('KriggingFinal.mat')
load('kriggedFinalPureOOS.mat')
for m=1:12
    mon = months{m};
    for h=0:23
        if m==9 && h<=11
            continue
        end
        f = sprintf('%sH%d',char(mon),h);
        ws = valsoos.(f);
        a1 = A(24*(m-1)+h+1);
        c1 = C(24*(m-1)+h+1);
        nugget = Nugget(24*(m-1)+h+1);
        
        fitsemi = @(h,a1,c1) sph(h,a1,c1); 
        fitcov = @(h,a1,c1,a2,c2) nugget-sph(h,a1,c1);
        krigged_val = zeros(size(ws,1),n);
        
        for d=1:size(ws,1)
            for p = 1:n
                %krigged_val(d,p) = kriging3(p,newlat, newlon, lats, lons, ws(d,:), fitcov, fitsemi, nugget, sph, a1, c1);
                krigged_val(d,p) = kriging(p,newlat, newlon, lats, lons, ws(d,:), fitcov, fitsemi, nugget, sph, '', a1, c1, 0, 0);
            end
        end
        kriggedFinal.(f) = krigged_val;
        sprintf('Done with Month %s, Hour %d',months{m},h)
    end
    save('kriggedFinalPureOOS.mat','kriggedFinal')
end

% save('KriggingFinal.mat','kriggedFinal')

% figure
% scatter(newlon,newlat, 20, krigged_alpha1, 'filled');
% colormap(hsv) %can be removed
% xlabel('Lon')
% ylabel('Lat')
% colorbar
% title('Kriging results for $\alpha_1$', 'Interpreter', 'latex')


%% Merge the files obtained
localFile = load('kriggedFinalPureOOS.mat');
terminalFile = load('terminal_kriggedFinalPureOOS.mat');


for m=1:12
    mon = months{m};
    for h=0:24
        f = sprintf('%sH%d',char(mon),h);
        isLocal = isfield(localFile.kriggedFinal,f);
        isTerminal = isfield(terminalFile.kriggedFinal,f);
        if isLocal
            kriggedFinal.(f) = localFile.kriggedFinal.(f);
        elseif isTerminal
            kriggedFinal.(f) = terminalFile.kriggedFinal.(f);
        end
    end
end
% 
save('kriggedFinal_Pureoos.mat','kriggedFinal')

%% Arrange the timeseries within
load('kriggedFinal_Pureoos.mat','kriggedFinal')
totalKrigged = nan(size(valuesoos,1),n);
for h=0:23
    idxJan = find(month(datesoos)==1 & hour(datesoos)==h);
    idxFeb = find(month(datesoos)==2 & hour(datesoos)==h);
    idxMar = find(month(datesoos)==3 & hour(datesoos)==h);
    idxApr = find(month(datesoos)==4 & hour(datesoos)==h);
    idxMay = find(month(datesoos)==5 & hour(datesoos)==h);
    idxJun = find(month(datesoos)==6 & hour(datesoos)==h);
    idxJul = find(month(datesoos)==7 & hour(datesoos)==h);
    idxAug = find(month(datesoos)==8 & hour(datesoos)==h);
    idxSep = find(month(datesoos)==9 & hour(datesoos)==h);
    idxOct = find(month(datesoos)==10 & hour(datesoos)==h);
    idxNov = find(month(datesoos)==11 & hour(datesoos)==h);
    idxDec = find(month(datesoos)==12 & hour(datesoos)==h);
    
    for m=1:12
        mon = months{m};
        f = sprintf('%sH%d',char(mon),h);
        idx = sprintf('idx%s',char(mon));
        totalKrigged(eval(idx),:) = kriggedFinal.(f);
    end
end
%% Split the total result in 5 years (for dimension handling)
idx = year(datesin)==2015;
datesin2015 = datesin(idx);
totalKrigged2015 = totalKrigged(idx,:);
idx = year(datesin)==2016;
datesin2016 = datesin(idx);
totalKrigged2016 = totalKrigged(idx,:);
idx = year(datesin)==2017;
datesin2017 = datesin(idx);
totalKrigged2017 = totalKrigged(idx,:);
idx = year(datesin)==2018;
datesin2018 = datesin(idx);
totalKrigged2018 = totalKrigged(idx,:);
idx = year(datesin)==2019;
datesin2019 = datesin(idx);
totalKrigged2019 = totalKrigged(idx,:);

%% Give the kriging maps the density we actually want

oldlat = newlat;
oldlon = newlon;

% The density we want
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


%  Intersect

idxNewtoOld = find(ismember([newlat newlon],[oldlat oldlon], 'rows'));

krigged_result2015 = nan(length(newlon),size(totalKrigged2015,1));
krigged_result2016 = nan(length(newlon),size(totalKrigged2016,1));
krigged_result2017 = nan(length(newlon),size(totalKrigged2017,1));
krigged_result2018 = nan(length(newlon),size(totalKrigged2018,1));
krigged_result2019 = nan(length(newlon),size(totalKrigged2019,1));
for i=1:size(totalKrigged2015,1)
    krigged_result2015(idxNewtoOld,i) = totalKrigged2015(i,:);
end
for i=1:size(totalKrigged2016,1)
    krigged_result2016(idxNewtoOld,i) = totalKrigged2016(i,:);
end
for i=1:size(totalKrigged2017,1)
    krigged_result2017(idxNewtoOld,i) = totalKrigged2017(i,:);
end
for i=1:size(totalKrigged2018,1)
    krigged_result2018(idxNewtoOld,i) = totalKrigged2018(i,:);
end
for i=1:size(totalKrigged2019,1)
    krigged_result2019(idxNewtoOld,i) = totalKrigged2019(i,:);
end
% Fill up the gaps

for i=1:size(totalKrigged2015,1)
    F = scatteredInterpolant(oldlon,oldlat,totalKrigged2015(i,:)','linear');
    krigged_result2015(:,i) = F(newlon, newlat);
end
for i=1:size(totalKrigged2016,1)
    F = scatteredInterpolant(oldlon,oldlat,totalKrigged2016(i,:)','linear');
    krigged_result2016(:,i) = F(newlon, newlat);
end
for i=1:size(totalKrigged2017,1)
    F = scatteredInterpolant(oldlon,oldlat,totalKrigged2017(i,:)','linear');
    krigged_result2017(:,i) = F(newlon, newlat);
end
for i=1:size(totalKrigged2018,1)
    F = scatteredInterpolant(oldlon,oldlat,totalKrigged2018(i,:)','linear');
    krigged_result2018(:,i) = F(newlon, newlat);
end
for i=1:size(totalKrigged2019,1)
    F = scatteredInterpolant(oldlon,oldlat,totalKrigged2019(i,:)','linear');
    krigged_result2019(:,i) = F(newlon, newlat);
end




% Make some plots
myidx = find(datesin2018==datenum('22-Dec-2018 02:00:00'));
figure
scatter(newlon,newlat, 20, krigged_result2018(:, myidx), 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results: 21-Dec-2018 18:00:00')

save('Results_Kriging_WindSpeeds.mat','datesin2015', 'krigged_result2015', 'datesin2016', 'krigged_result2016', ...
    'datesin2017', 'krigged_result2017', 'datesin2018', 'krigged_result2018', 'datesin2019', 'krigged_result2019', '-v7.3');
save('Results_Kriging_WindSpeeds2015.mat','datesin2015', 'krigged_result2015','-v7.3');
save('Results_Kriging_WindSpeeds2016.mat','datesin2016', 'krigged_result2016','-v7.3');
save('Results_Kriging_WindSpeeds2017.mat','datesin2017', 'krigged_result2017','-v7.3');
save('Results_Kriging_WindSpeeds2018.mat','datesin2018', 'krigged_result2018','-v7.3');
save('Results_Kriging_WindSpeeds2019.mat','datesin2019', 'krigged_result2019','-v7.3');
%% Make averages over a month
% March 2017
idxMar = month(datesin2017)==3;
valsMar = mean(krigged_result2017(:,idxMar),2);

figure
scatter(newlon,newlat, 20, valsMar, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results: March 2017')

% July 2017
idxJul = month(datesin2017)==7;
valsJul = mean(krigged_result2017(:,idxJul),2);

figure
scatter(newlon,newlat, 20, valsJul, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results: July 2017')

% September 2017
idxSep = month(datesin2017)==9;
valsSep = mean(krigged_result2017(:,idxSep),2);

figure
scatter(newlon,newlat, 20, valsSep, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results: September 2017')

% December 2017
idxDec = month(datesin2017)==12;
valsDec = mean(krigged_result2017(:,idxDec),2);

figure
scatter(newlon,newlat, 20, valsDec, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results: December 2017')



%Make maps for Pure OOS and OOS-PureOOS difference
pureoos = load('Results_Kriging_PureOOS_WindSpeeds.mat','krigged_result');
oos = load('Results_Kriging_OOS_WindSpeeds.mat','krigged_result');

% July 2019
idxJul = month(datesoos)==7;
valsJul = mean(krigged_result(:,idxJul),2);
valsJuldiff = mean(oos.krigged_result(:,idxJul)-pureoos.krigged_result(:,idxJul),2);

figure
scatter(newlon,newlat, 20, valsJul, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results derived out-of-sample: July 2019')

figure
scatter(newlon,newlat, 20, valsJuldiff, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Delta wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results out-of-sample difference: July 2019')





% September 2019
idxJul = month(datesoos)==9;
valsJul = mean(krigged_result(:,idxJul),2);
valsJuldiff = mean(oos.krigged_result(:,idxJul)-pureoos.krigged_result(:,idxJul),2);

figure
scatter(newlon,newlat, 20, valsJul, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results derived out-of-sample: September 2019')

figure
scatter(newlon,newlat, 20, valsJuldiff, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results out-of-sample difference: September 2019')


% December 2019
idxJul = month(datesoos)==12;
valsJul = mean(krigged_result(:,idxJul),2);
valsJuldiff = mean(oos.krigged_result(:,idxJul)-pureoos.krigged_result(:,idxJul),2);


figure
scatter(newlon,newlat, 20, valsJul, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results derived out-of-sample: December 2019')


figure
scatter(newlon,newlat, 20, valsJuldiff, 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
a=colorbar;
ylabel(a,'Wind speed (m/s)','FontSize',9,'Rotation',270);
a.Label.Position(1) = 3;
%colorbar
title('Kriging results out-of-sample difference: December 2019')


%% Function to determine the variogram and covariance functions
% function [g, cf] = calculate_semivariogram(z, edges, D, B, angle)
%     for k=1:length(edges)-1
%         [idxrow, idxcol] = find(D>=edges(k) & D<edges(k+1) & ((B>=angle-20 & B<=angle+20) | B==0));
% 
%         gamma = 0;
%         for i=1:length(idxrow)
%             gamma = gamma + (z(idxrow(i))-z(idxcol(i)))^2;
%         end
%         num = length(idxrow);
%         g(k) = gamma/(2*num);
%         C = cov(z(idxrow),z(idxcol));
%         cf(k) = C(1,2);
%     end
%     k = length(edges);
%     [idxrow, idxcol] = find(D>=edges(k-1) & ((B>=angle-20 & B<=angle+20) | B==0));
%     gamma = 0;
%     for i=1:length(idxrow)
%         gamma = gamma + (z(idxrow(i))-z(idxcol(i)))^2
%     end
%     C = cov(z(idxrow),z(idxcol));
%     cf(k) = C(1,2);
%     num = length(idxrow);
%     g(k) = gamma/(2*num);
% end
% 


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

