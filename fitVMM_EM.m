function obj =  fitVMM_EM(X,Kmax,varargin)
[start,tol,maxIter,Regularize,debg] = process_options(varargin,'start',[],'tol',1e-12,'maxIter',300,'Regularize',1e-3,'debg',2);
[N,D]= size(X);
%% Initialization using Spherical KMeans
% Inititalize mean, kappa and component mixing weight
if isempty(start)
    [Mu,Kappa,W] = InitParameters(X,Kmax,'kappaMin',Regularize,'debg',0);
else
    W = start.W;
    Mu = start.MU;
    Kappa = start.Kappa;
    Kmax = length(W);
end
% EM
iteration = 2;
loglike = -inf(maxIter,1);
converged = 0;
while (~converged)
    % E-step
    [R,loglike(iteration)] = Expectation(X,W,Mu,Kappa);
    % M-step
    [W,Mu,Kappa] = Maximization(X,Kappa,R);
    Kappa = Kappa + 1/Regularize;
    % convergence check
    deltlike = loglike(iteration) - loglike(iteration-1);
    deltlike = abs(100*(deltlike/loglike(iteration-1)));
    if(deltlike < tol || iteration > maxIter)
        converged = 1;
        loglike(iteration+1:end) = [];
    end
    prt(debg,1,sprintf('########### EM Iteration: %d, LogLikelihood=%8.8f, Delta=',iteration,loglike(iteration)),deltlike);
    iteration = iteration + 1;
end
% Store results in object
obj.Iters = iteration-1;
obj.DistName = 'Mixture of vMF distributions';
obj.NDimensions = D;
obj.NSamples = N;
obj.NComponents = Kmax;
obj.PComponents = W;
obj.mu = Mu;
obj.Kappa = Kappa;
obj.E = R;
[~, idx] = max(R,[],2);
obj.Class =  idx;
obj.logL =  loglike(end);
end

function [R,llh] =  Expectation(X,W,Mu,Kappa)
% E-Step
[N, D]= size(X);
logNormalize  = log(W) + (D/2-1)*log(Kappa) - (D/2)*log(2*pi) - logbesseli(D/2-1,Kappa);
R  = X * (Mu'.*(ones(D,1)*Kappa));
R  = bsxfun(@plus,R,logNormalize);
T = logsumexp(R,2);
llh = sum(T)/N; % loglikelihood
R = exp(bsxfun(@minus,R,T));
end

function [W,Mu,Kappa] = Maximization(X,Kappa,R)
[N, K] = size(R);
[~, D] = size(X);
W = sum(R,1)./N;
Mu = R'*X;
for k=1:K
    normMu   = sqrt(Mu(k,:)*Mu(k,:)');
    rbar  = normMu/W(k);
    Mu(k,:)  = Mu(k,:)/normMu;
    Kappa(k) = (rbar*D - rbar^3)/(1-rbar^2);
end
end

function [logb] = logbesseli(nu,x)
% log of the Bessel function, extended for large nu and x
% approximation from Eqn 9.7.7 of Abramowitz and Stegun
% http://www.math.sfu.ca/~cbm/aands/page_378.htm
frac = x/nu;
square = 1 + frac.^2;
root = sqrt(square);
eta = root + log(frac) - log(1+root);
approx = - log(sqrt(2*pi*nu)) + nu*eta - 0.25*log(square);
logb = approx;
% alternative less accurate approximation from Eqn 9.6.7
% this gives better results on the Classic400 dataset!
% logb = nu*log(x/2) - gammaln(nu+1);
% [bessel,flags] = besseli(nu,x);
% if any(flags ~= 0) || any(bessel == Inf)
%     besselproblem = [x, bessel, flags];
% end
% logb = bessel;
% nz = find(bessel > 0);
% z = find(bessel == 0);
% logb(nz) = log(bessel(nz));
% logb(z) = nu*log(x(z)/2) - gammaln(nu+1);
%[nu*ones(size(x))'; x'; approx'; logb']
end

function prt(debg, level, txt, num)
% Print text and number to screen if debug is enabled.
if(debg >= level)
    if(numel(num) == 1)
        disp([txt num2str(num)]);
    else
        disp(txt)
        disp(num)
    end
end
end
