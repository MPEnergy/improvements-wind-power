function krigged_alpha = kriging(p, newlat, newlon, lats, lons, alpha, covfun, fitalpha, nugget, fun1, fun2, param1, param2, param3, param4)
if  isempty(find(lats==newlat(p))) || isempty(find(lons(find(lats==newlat(p)))==newlon(p))) %    newlat(p)-fix(newlat(p))~=0.5 || newlon(p)-fix(newlon(p))~=0.5
    %C*lambda=v
    n = length(lats);
    v = [];
    for i=1:n
        [dist,~] = lldistkm([lats(i) lons(i)], [newlat(p) newlon(p)]);
        v = [v covfun(dist, param1, param2, param3, param4)];
    end
    
    for i=1:n
        for j=1:n
            [dist,~] = lldistkm([lats(i) lons(i)], [lats(j) lons(j)]);
            C(i,j) = covfun(dist, param1, param2, param3, param4);
            if i==j
                C(i,j) = nugget;%C(i,j) = 0.0033;
            end
        end
    end
    % Add the unbiased constraint: that all the weights sum up to 1
    C = [C ones(n,1); ones(1,n) 0];
    sol = inv(C)*[v';1];
    lambda = sol(1:end-1);
    
    aux = 0;
    for i=1:n
        aux = aux+lambda(i)*(alpha(i)-mean(alpha));
    end
    krigged_alpha = aux+ mean(alpha) ;
   
    
else %it is already a point in the grid
    idx = find(lats==newlat(p) & lons==newlon(p));
    
    if ~isempty(idx)
        krigged_alpha = alpha(idx);
    else
        sprintf('Wrong lat %d, lon %d', newlat(p), newlon(p));
        return;
    end
end
end