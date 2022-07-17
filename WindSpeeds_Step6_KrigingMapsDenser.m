%% Make the cokriging maps denser: for curve kriging

%% 1. Read the cokriging profiles
clear all
p = xlsread('SeasonalityResiduals.xlsx');

lats = p(:,2);
lons = p(:,3);
clear p
%load('CokrigingResults_SigmaBSB.mat')
load('CokrigingResults_SigmaBSB_NegCorrectedComplete.mat')
oldlat = newlat;
oldlon = newlon;

%% 2. The density we want
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

%% 3. Intersect

idxNewtoOld = find(ismember([newlat newlon],[oldlat oldlon], 'rows'));

krigged_result = nan(length(newlon),size(krigged_total,2));
for i=1:size(krigged_total,2)
    krigged_result(idxNewtoOld,i) = krigged_total(:,i);
end

%% 4. Fill up the gaps

for i=1:size(krigged_total,2)
    F = scatteredInterpolant(oldlon,oldlat,krigged_total(:,i),'linear');
    krigged_result(:,i) = F(newlon, newlat);
end


figure
scatter(newlon,newlat, 20, krigged_result(:,283), 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar
title('Kriging results for seasonality in December, hour 18.00')


figure
scatter(oldlon,oldlat, 20, krigged_total(:,1), 'filled');
colormap(hsv) %can be removed
xlabel('Lon')
ylabel('Lat')
colorbar

krigged_sigmaBSB = krigged_result;
save('CokrigingResults_SigmaBSBDense_NegativeValuesRemovedComplete.mat','newlat','newlon','krigged_sigmaBSB')
