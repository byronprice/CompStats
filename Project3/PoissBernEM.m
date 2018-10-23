function [piParam,lambda1,lambda2] = PoissBernEM(data)
%PoissBernEM.m, Project 3, 1-c
%   Fit latent variable Poisson-Bernoulli model for count data using the EM
%    algorithm. Data is assumed to have either a Poisson(lambda1) or a
%    Poisson(lambda2), depending on a Bernoulli distribution with parameter
%    piParam
%     e.g. data ~ Poisson(lambda1) with probability piParam and
%          data ~ Poisson(lambda2) with probability (1-piParam)
%INPUT: data - count data of observations, e.g. number of accidents at a
%         set of traffic intersections on a given day
%OUTPUTS:
%       piParam - probability that observed data point comes from Poisson
%          with parameter lambda1
%       lambda1 - Poisson distribution free parameter lambda for group 1
%       lambda2 - Poisson distribution free parameter lambda for group 2
%
%Created: 2018/10/22
% By: Byron Price

% ensure data is a column vector
data = data(:);

% initialize parameters
piParam = 0.5;
N = length(data);
dataMean = mean(data);
lambda1 = sum(data.*(data>dataMean))/sum(data>dataMean);
lambda2 = sum(data.*(data<=dataMean))/sum(data<=dataMean);

maxIter = 1e5;
tolerance = 1e-8;
for tt=1:maxIter
    % E step, calculate alpha-i,t
    tmp = poisspdf(data,lambda1)*piParam;
    tmp2 = poisspdf(data,lambda2)*(1-piParam);
    alpha = tmp./(tmp+tmp2);
    
    % M step, calculate new values for parameters pi, lambda1, lambda2
    alphaSum = sum(alpha);
    piParamStar = alphaSum/N;
    lambda1Star = sum(alpha.*data)/alphaSum;
    lambda2Star = sum((1-alpha).*data)/sum(1-alpha);
    
    totalPrecision = abs(piParamStar-piParam)+abs(lambda1Star-lambda1)+abs(lambda2Star-lambda2);
    if totalPrecision<=tolerance
        break;
    end
    
    piParam = piParamStar;
    lambda1 = lambda1Star;
    lambda2 = lambda2Star;
end
end

