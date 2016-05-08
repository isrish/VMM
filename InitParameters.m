function [mu,kappa,w] = InitParameters(X,K,varargin)
%% InitParameters
% Initialize Mean Vector and Concentration parameters of Mixture of von
% Mises-Fisher distributions. First the points in the N-by-D data matrix X are partitioned into K clusters
% using Spherical Kmeans.
% 
%   [mu,kappa,w] = InitParameters(X, K) returns the cluster centroid locations in the mu (K-by-D matrix center.
%   [mu,kappa,w] =  InitParameters(..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%   optional parameter name/value pairs to control the iterative algorithm used by KMEANS.  
%   Parameters are:
%
%   'kappaMin'  -  the minimum concentration parameter; Default is 1000.
%   'reps' - Number of times to repeat the clustering, each with a new set of initial centroids. Default is 10. 
%   'debg' = Verbosity level.  0 [default]

%#   $Author: Israel D. Gebru $    $Date: 2016/04/28 $    $Revision: 1.0 $
%#   Copyright:

%% Initialize parameters using Spherical Kmeans 
[kappaMin, debg,reps] = process_options(varargin,'kappaMin',100,'debg',0,'reps',10);
[n,d] = size(X);
%%
bestlabel = [];
sumD = zeros(1,K);
bCon = false;
maxit = 100;
%%
for t=1:reps
    mu = X(randsample(n,K),:);
    last = 0;label=1;
    it = 0;
    while any(label ~= last) && it<maxit
        last = label;
        prt(debg,2,sprintf('########### Initialization Iteration: '),it);
        % assign points to nearest cluster
        simMat =  X*mu';
        [val,label] =  max(simMat,[],2);
        ll = unique(label);
        if length(ll) < K
            missCluster = 1:K;
            missCluster(ll) = [];
            missNum = length(missCluster);
            [~,idx] = sort(val);
            label(idx(1:missNum)) = missCluster;
        end
        E = sparse(1:n,label,1,n,K,n);  % transform label into indicator matrix
        mu = full((E*spdiags(1./sum(E,1)',0,K,K))'*X);    % compute center of each cluster
        centernorm = sqrt(sum(mu.^2, 2));
        mu = mu ./ centernorm(:,ones(1,d));
        it=it+1;
    end
    if it<maxit
        bCon = true;
    end
    if isempty(bestlabel)
        bestlabel = label;
        bestcenter = mu;
        if reps>1
            if any(label ~= last)
                simMat = full(X*mu');
            end
            D = 1-simMat;
            for k = 1:K
                sumD(k) = sum(D(label==k,k));
            end
            bestsumD = sumD;
            bestD = D;
        end
    else
        if any(label ~= last)
            simMat=full(X*mu');
        end
        D = 1-simMat;
        for k = 1:K
            sumD(k) = sum(D(label==k,k));
        end
        if sum(sumD) < sum(bestsumD)
            bestlabel = label;
            bestcenter = mu;
            bestsumD = sumD;
            bestD = D;
        end
    end
end
% initializing kappa, mixing weight 
mu = bestcenter;
kappa = kappaMin*ones(1,K);
w = ones(1,K).*(1/n);
for k = 1:K
    idx = bestlabel==k;
    if sum(idx)>0
        w(k) = sum(idx)/n;        
        normMu   = sqrt(mu(k,:)*mu(k,:)');
        rbar  = normMu/w(k);
        mu(k,:)  = mu(k,:)./normMu;
        kappa(k) = max((rbar*d - rbar^3)/(1-rbar^2), kappaMin);
    else
        prt(debg,3,'Empty Cluster found!',[]);
    end
end
kappa = kappa * d;
w = w./sum(w);
end