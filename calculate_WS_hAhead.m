function newWS = calculate_WS_hAhead(thisLat, thisLon, allLat, allLon, params, dates, ws, withError, hoursAhead, randGenerator)
%% Function to simulate N paths of the wind speed at a new location and 95%PI

%% Find the closest (up to 2 digits) lat and lon
idx = find(allLat==round(thisLat*100/100,2) & allLon==round(thisLon*100/100,2));
m = params.m(idx);
alpha1 = params.alpha1(idx);    alpha2 = params.alpha2(idx);
beta1 = params.beta1(idx);      beta2 = params.beta2(idx);      beta3 = params.beta3(idx);      beta4 = params.beta4(idx);
beta5 = params.beta5(idx);      beta6 = params.beta6(idx);      beta7 = params.beta7(idx);      beta8 = params.beta8(idx);
beta9 = params.beta9(idx);      beta10 = params.beta10(idx);    beta11 = params.beta11(idx);    beta12 = params.beta12(idx);
beta13 = params.beta13(idx);    beta14 = params.beta14(idx);    beta15 = params.beta15(idx);    beta16 = params.beta16(idx);
beta17 = params.beta17(idx);    beta18 = params.beta18(idx);    beta19 = params.beta19(idx);    beta20 = params.beta20(idx);
beta21 = params.beta21(idx);    beta22 = params.beta22(idx);    beta23 = params.beta23(idx);    beta24 = params.beta24(idx);
betas = [beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10,beta11,beta12,beta13,beta14,beta15,beta16,beta17,beta18,beta19,beta20,beta21,beta22,beta23,beta24];
a0 = params.a0(idx);    a1 = params.a1(idx);    a2 = params.a2(idx);
a3 = params.a3(idx);    a4 = params.a4(idx);    a5 = params.a5(idx);
a6 = params.a6(idx);    a7 = params.a7(idx);    a8 = params.a8(idx);
a9 = params.a9(idx);    a10 = params.a10(idx);  a11 = params.a11(idx);
a12 = params.a12(idx);
p = params.p(idx,:)';
q = params.q(idx,:)';
delta = params.delta(idx);    gamma0 = params.gamma0(idx);    gamma1 = params.gamma1(idx);    lambda = params.lambda(idx);
omega = params.omega(idx);    skew = params.skew(idx);        shape = params.shape(idx);
%% Check invertibility of the ARMA(2,24) model, make the model invertible if it is not
if hoursAhead<24
    %[betas, noiseVarianceMultiplier] = local_checkInvertibilityARMA(alpha1, alpha2, betas(hoursAhead:24));
        %betas = [zeros(1, hoursAhead-1) betas(hoursAhead:24)];
        [betas, noiseVarianceMultiplier] = local_checkInvertibilityARMA(alpha1, alpha2, betas);

else
    [betas, noiseVarianceMultiplier] = local_checkInvertibilityARMA(alpha1, alpha2, betas);
   
end

if noiseVarianceMultiplier~=1
    betas = fliplr(betas);
end

% if noiseVarianceMultiplier~=1
%     beta24 = betas(1); beta23 = betas(2); beta22 = betas(3); beta21 = betas(4); beta20 = betas(5); beta19 = betas(6);
%     beta18 = betas(7); beta17 = betas(8); beta16 = betas(9); beta15 = betas(10); beta14 = betas(11); beta13 = betas(12);
%     beta12 = betas(13); beta11 = betas(14); beta10  = betas(15); beta9 = betas(16); beta8 = betas(17); beta7 = betas(18);
%     beta6 = betas(19); beta5 = betas(20); beta4 = betas(21); beta3 = betas(22); beta2 = betas(23); beta1 = betas(24);
%     
%     betas = [beta1,beta2,beta3,beta4,beta5,beta6,beta7,beta8,beta9,beta10,beta11,beta12,beta13,beta14,beta15,beta16,beta17,beta18,beta19,beta20,beta21,beta22,beta23,beta24];
% end
% Persistence of the model should be <1
if gamma0+gamma1/2+delta>=1
    return
end
%% Start the simulation
ws = ws.^m; %comment this out when doing time validation

%sigma_GARCH
sigma_GARCH = nan(length(dates),1);
sigma_GARCH(1) = nanstd(ws); %
sigma_GARCH(1) = omega/(1-gamma0-gamma1/2-delta);
%sigma_BSB
IdxVec = (month(dates)-ones(size(dates))).*24+hour(dates)+ones(size(dates));
sigma_BSB = sqrt(q(IdxVec));

%Xi
if ~withError
    Xi = zeros(size(dates));
else
    if isempty(randGenerator)
        randGenerator = rng;
    end
    rng(randGenerator);
    Xi = skewtdis_rnd(shape,skew,length(dates));
end





%% S: seasonality part in the mean
vv = round(dates.*24)./24-datenum(year(dates),1,1);
S = repmat(a0,size(dates));
for i =1:6
    varcos = eval(sprintf('a%d',2*i-1));
    varsin = eval(sprintf('a%d',2*i));
    S = S+ varcos*cos(2*i*pi*vv/365.25)+varsin*sin(2*i*pi*vv/365.25);
end
S = S+p(IdxVec);

%S = hpfilter(S,1);


%% Start the simulation
%epsilon
epsilon = nan(length(dates),1);
% epsilon(1:hoursAhead) = 0;
% epsilontilde(1:hoursAhead) = epsilon(1:hoursAhead);
% WS(1:hoursAhead) = ws(1:hoursAhead);

% Test
epsilon(1:24) = 0;
epsilontilde(1:24) = epsilon(1:24);
WS(1:24) = ws(1:24);
for t=2:24
    sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
    
    epsilon(t) = sigma_GARCH(t)*Xi(t)*noiseVarianceMultiplier;
    
    epsilontilde(t) = lambda*epsilontilde(t-1)+epsilon(t);
    
    epsilonhat(t) = sqrt(sigma_BSB(t))*epsilontilde(t);
    
end


% for t=2:hoursAhead
%     sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
%     
%     epsilon(t) = sigma_GARCH(t)*Xi(t)*noiseVarianceMultiplier;
%     
%     epsilontilde(t) = lambda*epsilontilde(t-1)+epsilon(t);
%     
%     epsilonhat(t) = sqrt(sigma_BSB(t))*epsilontilde(t);
%     
% end



% for t = hoursAhead+1:24
%     sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
%     epsilon(t) = sigma_GARCH(t)*Xi(t)*noiseVarianceMultiplier;
%     epsilontilde(t) = lambda*epsilontilde(t-1)+epsilon(t);
%     epsilonhat(t) = sqrt(sigma_BSB(t))*epsilontilde(t);
%     if hoursAhead>2
%         WS(t) = S(t) + alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2));%+epsilonhat(t);
%         for k=hoursAhead:24
%             if t-k>0
%                 WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
%             else
%                 continue
%             end
%             
%         end
% %         k = 1;
% %         while t-hoursAhead>0 && k>=hoursAhead
% %             WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
% %             k = k-1;
% %         end
%     elseif hoursAhead==2
%         WS(t) = S(t) + alpha1*(WS(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+epsilonhat(t)+betas(1)*epsilonhat(t-1);
%     else
%         WS(t) = S(t) + alpha1*(wS(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+epsilonhat(t);
%     end
%     
% epsilonhat(t) = ws(t)-WS(t);%-ws(t);
% epsilontilde(t) = epsilonhat(t)/sqrt(sigma_BSB(t));
% epsilon(t) = epsilontilde(t)-lambda*epsilontilde(t-1);
% sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
% end

% epsilonhat(24) = WS(24)-ws(24);
% epsilontilde(24) = epsilonhat(24)/sqrt(sigma_BSB(24));
% epsilon(24) = epsilontilde(24)-lambda*epsilontilde(23);
% sigma_GARCH(24) = sqrt(omega+(gamma0+gamma1*(epsilon(23)>0))*epsilon(23)^2+delta*sigma_GARCH(23)^2);

% % % % % % for t=25:length(dates)
% % % % % %     sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
% % % % % %     epsilon(t) = sigma_GARCH(t)*Xi(t)*noiseVarianceMultiplier;
% % % % % %     epsilontilde(t) = lambda*epsilontilde(t-1)+epsilon(t);
% % % % % %     epsilonhat(t) = sqrt(sigma_BSB(t))*epsilontilde(t);
% % % % % %     
% % % % % %    
% % % % % %     
% % % % % %     if hoursAhead>2
% % % % % %         WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2));%+epsilonhat(t)*(withError);
% % % % % %         if hoursAhead<=24
% % % % % %             for k=hoursAhead:24
% % % % % %                 if t-k>0
% % % % % %                     WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
% % % % % %                 else
% % % % % %                     continue
% % % % % %                 end
% % % % % %                 
% % % % % %             end
% % % % % %         end
% % % % % %     elseif hoursAhead==2
% % % % % %         WS(t) = S(t) + alpha1*(WS(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+epsilonhat(t)+betas(1)*epsilonhat(t-1);
% % % % % %     else
% % % % % %         WS(t) = S(t) + alpha1*(ws(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+epsilonhat(t);
% % % % % %     end
% % % % % % %     WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2))+beta1*epsilonhat(t-1)+beta2*epsilonhat(t-2)+...
% % % % % % %             beta3*epsilonhat(t-3)+beta4*epsilonhat(t-4)+beta5*epsilonhat(t-5)+beta6*epsilonhat(t-6)+beta7*epsilonhat(t-7)+...
% % % % % % %             beta8*epsilonhat(t-7)+beta9*epsilonhat(t-9)+beta10*epsilonhat(t-10)+beta11*epsilonhat(t-11)+beta12*epsilonhat(t-12)+...
% % % % % % %             beta13*epsilonhat(t-13)+beta14*epsilonhat(t-14)+beta15*epsilonhat(t-15)+beta16*epsilonhat(t-16)+beta17*epsilonhat(t-17)+...
% % % % % % %             beta18*epsilonhat(t-18)+beta19*epsilonhat(t-19)+beta20*epsilonhat(t-20)+beta21*epsilonhat(t-21)+beta22*epsilonhat(t-22)+...
% % % % % % %             beta23*epsilonhat(t-23)+beta24*epsilonhat(t-24)+epsilonhat(t)*(withError);     
% % % % % %     epsilonhat(t) = ws(t)-WS(t);
% % % % % %     epsilontilde(t) = epsilonhat(t)/sqrt(sigma_BSB(t));
% % % % % %     epsilon(t) = epsilontilde(t)-lambda*epsilontilde(t-1);
% % % % % %     sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
% % % % % % end

%Test
for t=25:length(dates)
    sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
    epsilon(t) = sigma_GARCH(t)*Xi(t)*noiseVarianceMultiplier;
    epsilontilde(t) = lambda*epsilontilde(t-1)+epsilon(t);
    epsilonhat(t) = sqrt(sigma_BSB(t))*epsilontilde(t);
    
    if t>=hoursAhead
        WS(t) = local_calculatePrediction(t, hoursAhead, alpha1, alpha2, betas, lambda, WS, ws, S, epsilonhat, epsilontilde, sigma_BSB);
    else
        WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2));
    end
    
    
    epsilonhat(t) = ws(t)-WS(t);
    epsilontilde(t) = epsilonhat(t)/sqrt(sigma_BSB(t));
    epsilon(t) = epsilontilde(t)-lambda*epsilontilde(t-1);
    sigma_GARCH(t) = sqrt(omega+(gamma0+gamma1*(epsilon(t-1)>0))*epsilon(t-1)^2+delta*sigma_GARCH(t-1)^2);
end
WS(WS<0) = 0;
WS(WS>(25^m)) = 25^m;


newWS = nthroot(WS, m)';
%%










%sigma = sigma_GARCH.^2.*sigma_BSB.*noiseVarianceMultiplier^2;
%[epsilon,sigma_GARCH] = local_calculateEpsilon(epsilon,alpha1,alpha2,betas,delta0,gamma1,gamma2,S,sigma_GARCH,sigma_BSB,ws,t);

%% Determine the Var(epsilon_t)=sigma^2'sum(Psi^2), Psi, coefficients from
%the Ma(infty) development of the ARMA(2,24) model

%Psi = local_MArepresentation(alpha1, alpha2, betas);


%sigma_epsilon = sigma*sum(Psi.^2);



% ws0 = repmat(ws(1),1, N);
% ws1 = repmat(ws(2),1, N);
% epsilon(1) = 0;
% epsilon(2) = 0;
% % vv = round(dates.*24)./24-datenum(year(dates),1,1);
% vv = DayOfYear(dates);
% if isempty(randGenerator)
%     randGenerator = rng;
% end
% rng(randGenerator);
% Xi = randn(length(dates),N);
% Xi = (Xi-mean(Xi))./std(Xi);
% for d = 3:length(dates)
%     %(1) Seasonality term
%     S = local_seasonality(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,vv(d));
%     s1 = local_seasonality(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,vv(d-1));
%     s0 = local_seasonality(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,vv(d-2));
%     %(2) AR(2) terms
%     A = alpha1*(ws1-s1)+alpha2*(ws0-s0);
%     
%     %(3) Residual terms
%     sigma = local_sigma(b0,b1,b2,vv(d));
%     %xi = randn(N,1);
%     if withErr
% %         xi = randn(N,1);
%         xi = Xi(d);
%     else
%         xi = 0;
%     end
%     
%     epsilon(d) = sqrt(sigma).*xi;
%     
%     newWS(d,:) = (S+A'+epsilon(d))'; 
%     ws0 = ws1;
% %     ws1 = newWS(d,:);
%     ws1 = ws(d); 
% end



end

function [epsilon,sigma_GARCH] = local_calculateEpsilon(epsilon,alpha1,alpha2,betas,delta0,gamma1,gamma2,S,sigma_GARCH,sigma_BSB,ws,tmaxknown)

%Assume the first 24 W's are known, calculate epsilon
epsilon(1) = ws(1)-S(1);
epsilon(2) = ws(2)-S(2)-alpha1*(ws(1)-S(1))-betas(1)*epsilon(1);
epsilonhat = epsilon./sqrt(sigma_BSB);
sigma_GARCH(2) = sqrt(delta0+gamma1*epsilonhat(1)^2+gamma2*sigma_GARCH(1)^2);
for t=3:tmaxknown
    k = 1;
    epsilon(t) = ws(t)-S(t)-alpha1*(ws(t-1)-S(t-1))-alpha2*(ws(t-2)-S(t-2));
    while t-k>0 && k<=24
        epsilon(t) = epsilon(t)-betas(k)*epsilon(t-k);
        k = k+1;
    end
    epsilonhat = epsilon./sqrt(sigma_BSB);
    sigma_GARCH(t) = sqrt(delta0+gamma1*epsilonhat(t-1)^2+gamma2*sigma_GARCH(t-1)^2);
end

end


function [newbetas, varMultiplier] = local_checkInvertibilityARMA(alpha1, alpha2, betas)
% https://stats.stackexchange.com/questions/406204/arima-simulation-followed-by-modeling-in-r-produces-bad-estimates/432342#432342
%Function to check invertibility of the ARMA process
varMultiplier = 1;
R_MA = roots([fliplr(betas) 1]);
R_AR = roots([-alpha2 -alpha1 1]);

if any(abs(R_MA)==1) | any(abs(R_AR)<=1) | ~isempty(intersect(R_AR, R_MA))
    return
end

idx = abs(R_MA)<1;
if sum(idx)~=0
    % Invert those roots with absolute value inside the unit circle
    newR = [R_MA(~idx); 1./R_MA(idx)];
    
    newbetas = poly(newR);
    newbetas = newbetas(1:24)./newbetas(25);
    
    %newbetas = newbetas./prod(1./R_MA(idx));
    %newbetas = newbetas./prod(R_MA(idx));
    varMultiplier = real(prod(R_MA(idx)));
else
    newbetas = betas;
end
end

function newval = local_calculatePrediction(t, hoursAhead, alpha1, alpha2, betas, lambda, WS, ws, S, epsilonhat, epsilontilde, sigma_BSB)

if hoursAhead>2 && hoursAhead<24
    %hoursAhead predictions
    WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2))+sigma_BSB(t)*lambda^hoursAhead*epsilontilde(t-hoursAhead);
    for k=1:hoursAhead-1
        WS(t) = WS(t)+betas(k)*lambda^(hoursAhead-k)*epsilontilde(t-hoursAhead);
    end
    for k=hoursAhead:24
        if t-k>0
            WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
        else
            continue
        end
    end
elseif hoursAhead>24
    if t~=hoursAhead
        WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2))+sigma_BSB(t)*lambda^hoursAhead*epsilontilde(t-hoursAhead);
        k = 1;
        while t-k>0 && k<=24
            WS(t) = WS(t)+betas(k)*sigma_BSB(t-k)*lambda^(hoursAhead-k)*epsilontilde(t-hoursAhead);
            k = k+1;
        end
    else
        WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2));
    end
elseif hoursAhead==24
    %24-hours ahead predictions
    WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(WS(t-2)-S(t-2))+sigma_BSB(t)*lambda^24*epsilontilde(t-24);
    for k=1:23
        if t-k>0
            WS(t) = WS(t)+betas(k)*sigma_BSB(t-k)*lambda^(hoursAhead-k)*epsilontilde(t-hoursAhead);
        else
            continue
        end
    end
    if t>24
        WS(t) = WS(t)+betas(24)*epsilonhat(t-24);
    end
elseif hoursAhead==2
    %2-hours ahead predictions
    WS(t) = S(t)+alpha1*(WS(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+sigma_BSB(t)*lambda^2*epsilontilde(t-2)+betas(1)*lambda*epsilontilde(t-2);
    k = 2;
    while t-k>0 && k<=24
        WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
        k = k+1;
    end
else
    %1-hour ahead predictions
    WS(t) = S(t)+alpha1*(ws(t-1)-S(t-1))+alpha2*(ws(t-2)-S(t-2))+sigma_BSB(t)*lambda*epsilontilde(t-1);
    k = 1;
    while t-k>0 && k<=24
        WS(t) = WS(t)+betas(k)*epsilonhat(t-k);
        k = k+1;
    end
end
newval = WS(t);
end






function Psi = local_MArepresentation(alpha1, alpha2, betas)
%Solve system of linear equations for Psi

M = diag(ones(25,1))+diag(-alpha1*ones(24,1),-1)+diag(-alpha2*ones(23,1),-2);
b = [1 betas];


Psi = inv(M)*b';
end
