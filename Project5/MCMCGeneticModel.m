function [beta0,beta,theta,posterior] = MCMCGeneticModel(Y,X,G)
%MCMCGeneticModel.m   Project 5, 2
%   Fit model described in Problem 2 of Project 5 with hybrid Gibbs sampler 
%   MCMC. It's a model used for analysis of genetic data, which gives disease
%   status of individuals as a logistic regression. The goal is to use the
%   disease status and some information about genetic variability to infer 
%   which genes are relevant to the disorder.
%   
%INPUTS:
%        Y - disease state data
%        X - design matrix of genomic markers
%        G - graph of dependency relationships between markers
%OUTPUTS:
%        beta0 - offset for logistic regression
%        beta - regression coefficients for each marker
%        theta - indicator for association of marker with disease
%        posterior - log of posterior (up to proportionality constant)
% Byron Price, 2018/11/25

global tau0;
tau0 = 1e-4;
global tau;
tau = 1;
global h;
h = -100;
global Tinv;
Tinv = 1/100;
global neighborhoods;
global neighbors;

p = size(X,2);
N = size(Y,1);

neighborhoods = cell(p,2);
for ii=1:p
    tmp = [];
    for jj=1:p
        if G(ii,jj)==1
            tmp = [tmp;jj];
        end
    end
    neighborhoods{ii,1} = tmp;
    neighborhoods{ii,2} = length(tmp);
end

neighbors = zeros(1,2);
count = 0;
for ii=1:p
    for jj=1:ii-1
        if G(ii,jj)==1
            count = count+1;
            neighbors(count,1) = ii;
            neighbors(count,2) = jj;
        end
    end
end
global numNeighbors;
numNeighbors = count;clear count;

burnIn = 1e4;
numIter = 1e4+burnIn;
beta0 = zeros(1,numIter);beta = zeros(p,numIter);theta = zeros(p,numIter);

% use glm to initialize chain
[b,~,stats] = glmfit(X,Y,'binomial');
beta0(1) = b(1);beta(:,1) = b(2:end);theta(:,1) = stats.p(2:end)<0.2;
onesVec = ones(N,1);

posterior = zeros(numIter,1);
posterior(1) = LogPosterior(Y,X,beta0(1),beta(:,1),theta(:,1));
% figure;plot(1,posterior(1),'.');pause(1e-3);hold on;
for tt=2:numIter
    % gibbs step for each theta
    for jj=1:p
        tmpTheta = [theta(1:jj-1,tt);theta(jj:end,tt-1)];
        tmpTheta(jj) = 0;
        logprob0 = LogThetaConditional(tmpTheta,beta(jj,tt-1),jj);
        tmpTheta(jj) = 1;
        logprob1 = LogThetaConditional(tmpTheta,beta(jj,tt-1),jj);
        
        prob1 = exp(logprob1-LogSum([logprob0,logprob1],2));
        
        theta(jj,tt) = binornd(1,prob1);
    end
    
    % laplace-mh step for beta0
    [beta0(tt)] = LaplaceStepMH(beta0(tt-1),Y,onesVec,X*beta(:,tt-1),0,0,N);
    
    % laplace-mh step for beta
    for jj=1:p
        omega = 1./SigSquare(theta(jj,tt));
        tmpBeta = [beta(1:jj-1,tt);beta(jj:end,tt-1)];
        offset = beta0(tt)+X(:,1:jj-1)*tmpBeta(1:jj-1)+X(:,jj+1:end)*tmpBeta(jj+1:end);
        [beta(jj,tt)] = LaplaceStepMH(tmpBeta(jj),Y,X(:,jj),...
            offset,0,omega,N);
    end
    
    posterior(tt) = LogPosterior(Y,X,beta0(tt),beta(:,tt),theta(:,tt));
%     plot(tt,posterior,'.');pause(1e-3);
end

% eliminate samples from burn-in period
beta0 = beta0(burnIn+1:end);
beta = beta(:,burnIn+1:end);
theta = theta(:,burnIn+1:end);
posterior = posterior(burnIn+1:end);
end

function [posterior] = LogPosterior(Y,X,beta0,beta,theta)
mu = InvLogit(beta0+X*beta);
sigmasquare = SigSquare(theta);

posterior = sum(LogBernPDF(Y,mu))+sum(LogNormalPDF(beta,0,sigmasquare))+...
    LogThetaPDF(theta);
end

function [beta] = LaplaceStepMH(prevBeta,Y,x,offset,beta0,omega,N)
mu = InvLogit(x*prevBeta+offset);
S = sqrt(1./(sum(mu.*(1-mu).*x.^2)+omega));

starBeta = normrnd(prevBeta,S);
starMu = InvLogit(x*starBeta+offset);
starS = sqrt(1./(sum(starMu.*(1-starMu).*x.^2)+omega));

logPostConditional = LaplacePosterior(prevBeta,Y,x,offset,beta0,omega,N);
logPostConditionalStar = LaplacePosterior(starBeta,Y,x,offset,beta0,omega,N);

logR = logPostConditionalStar+LogNormalPDF(prevBeta,starBeta,starS.^2)-...
    (logPostConditional+LogNormalPDF(starBeta,prevBeta,S.^2));


if logR>=0 || log(rand)<logR
    beta = starBeta;
else
    beta = prevBeta;
end
end

function [logprob] = LaplacePosterior(beta,y,x,offset,beta0,omega,N)
eta = offset+x*beta;
for ii=1:N
    eta(ii) = y(ii).*eta(ii)-LogSum([0,eta(ii)],2);
end
logprob = sum(eta)-0.5*omega*(beta-beta0).^2;
end

% function [beta] = LaplaceStepMH_Beta(prevBeta,Y,x,offset,omega,X,beta0,index,theta)
% mu = InvLogit(x*prevBeta(index)+offset);
% S = sqrt(1./(sum(mu.*(1-mu).*x.^2)+omega));
% 
% starBeta = normrnd(prevBeta(index),S);
% % starMu = InvLogit(x*starBeta+offset);
% % starS = sqrt(1./(sum(starMu.*(1-starMu).*x.^2)+omega));
% 
% sigmasquare = SigSquare(theta(index));
% logPostConditional = sum(LogBernPDF(Y,InvLogit(beta0+X*prevBeta)))+...
%     LogNormalPDF(prevBeta(index),0,sigmasquare);
% newBeta = prevBeta;newBeta(index) = starBeta;
% logPostConditionalStar = sum(LogBernPDF(Y,InvLogit(beta0+X*newBeta)))+...
%     LogNormalPDF(starBeta,0,sigmasquare);
% 
% % logR = logPostConditionalStar+LogNormalPDF(prevBeta(index),starBeta,starS.^2)-...
% %     (logPostConditional+LogNormalPDF(starBeta,prevBeta(index),S.^2));
% 
% logR = logPostConditionalStar-logPostConditional;
% 
% if logR>=0 || log(rand)<logR
%     beta = starBeta;
% else
%     beta = prevBeta(index);
% end
% end

function [logprob] = LogNormalPDF(x,mu,sigmasquare)
logprob = -0.5.*log(sigmasquare)-(x-mu).*(x-mu)./(2.*sigmasquare);
end

function [logprob] = LogBernPDF(y,mu)
logprob = y.*log(mu)+(1-y).*log(1-mu);
end

function [logprob] = LogThetaPDF(theta)
global h;global Tinv;global numNeighbors;global neighbors;
logprob = (h/2)*sum(2.*theta-1);
for ii=1:numNeighbors
    logprob = logprob+Tinv*(2*theta(neighbors(ii,1))-1)*(2*theta(neighbors(ii,2))-1);
end
end

function [logprob] = LogThetaConditional(theta,betaj,index)
global neighborhoods;global Tinv;global h;
sigmasquare = SigSquare(theta(index));

logprob = LogNormalPDF(betaj,0,sigmasquare);
if theta(index)==1
    neighborhood = neighborhoods{index,1};
    neighbors = neighborhoods{index,2};
    
    summation = 0;
    for ii=1:neighbors
        summation = summation+(2*theta(neighborhood(ii))-1);
    end
    summation = 2*Tinv*summation;
    logprob = logprob+h+summation;
end
end

function [y] = InvLogit(x)
y = 1./(1+exp(-x));
end

function [y] = Logit(p)
y = log(p)-log(1-p);
end

function [sigmasquare] = SigSquare(theta)
global tau0;global tau;
sigmasquare = (1-theta).*tau0+theta.*tau;
end