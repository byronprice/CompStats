function [mu1,mu2,sigsquare,Z,piParam] = GaussChangepointEM(data)
%GaussChangepointEM.m, Project 3, 2-d
%   Fit latent variable changepoint Gaussian mixture model using EM
%    algorithm. Data is assumed to have either a Normal(mu1,sigsquare) or a
%    Normal(mu2,sigsquare), depending on whether it occurs before or after
%    the changepoint Z
%     e.g. data-i ~ Normal(mu1,sigsquare) if i<=Z
%          data-i ~ Normal(mu2,sigsquare) if i>Z
%INPUT: data - observed data of size N-by-1, following changepoint normal mixture
%         distribution
%OUTPUTS:
%       mu1 - normal distribution mean parameter for data occurring before
%        the changepoint
%       mu2 - normal distribution mean parameter for data occurring after
%        the changepoint
%       sigsquare - variance of the data
%       Z - best estimate of the changepoint, provided as a number from 1
%        to N-1
%       piParam - data used to estimate changepoint
%
%Created: 2018/10/22
% By: Byron Price

% ensure data is a column vector
data = data(:);

% initialize parameters
N = length(data);
logData = log(data);

mu1 = mean(data(1:floor(N/2)));
mu2 = mean(data(ceil(N/2):N));
sigsquare = (var(data(1:floor(N/2)))+var(data(ceil(N/2):N)))/2;

maxIter = 1e5;
tolerance = 1e-6;
for tt=1:maxIter
    % E step, calculate pi-i,t
    tmp = zeros(N,1); 
    for ii=1:N
       meanVec = [mu1*ones(ii,1);mu2*ones(N-ii,1)];
       tmp(ii) = -(data-meanVec)'*(data-meanVec)./(2*sigsquare);
    end
    
    logpi = tmp-LogSum(tmp,N);
    
    expectationLessThan = zeros(N,1); % get expectations, stored as log(expectation)
    expectationMoreThan = zeros(N,1);    
    for ii=1:N
        expectationLessThan(ii) = LogSum(logpi(ii:N),N-ii+1);
        expectationMoreThan(ii) = LogSum(logpi(1:ii-1),ii-1);
    end
    
    % M step, calculate new values for parameters mu1, mu2, sigsquare
    
    % uses a trick, the sum over i of x-i times y-i is equal to 
    %  exp{log[sum of exp{log(x-i)+log(y-i)}]}
    %  which is useful if log(x-i) are very large negative numbers
    %  because all x-i will be effectively zero
    mu1Star = exp(LogSum(expectationLessThan+logData,N)-LogSum(expectationLessThan,N));
    mu2Star = exp(LogSum(expectationMoreThan+logData,N)-LogSum(expectationMoreThan,N));
    
    tmp = zeros(N,1);
    for ii=1:N
       meanVec = [mu1Star*ones(ii,1);mu2Star*ones(N-ii,1)];
       tmp(ii) = log((data-meanVec)'*(data-meanVec));
    end
    
    sigsquareStar = (1/N)*exp(LogSum(logpi+tmp,N));
    
    totalPrecision = abs(mu1-mu1Star)+abs(mu2-mu2Star)+abs(sigsquareStar-sigsquare);
    if totalPrecision<=tolerance
        [~,Z] = max(logpi);
        piParam = exp(logpi);
        break;
    end
    
    mu1 = mu1Star;
    mu2 = mu2Star;
    sigsquare = sigsquareStar;
end
end

function [summation] = LogSum(vector,vectorLen)
if vectorLen==0
    summation = -Inf;
elseif vectorLen==1
    summation = vector(1);
else
    vector = sort(vector);
    summation = LogSumExpTwo(vector(1),vector(2));
    for ii=2:vectorLen-1
        summation = LogSumExpTwo(summation,vector(ii+1));
    end
end

end

function [y] = LogSumExpTwo(x1,x2)
check = x1>=x2;
if check==1
    y = x1+SoftPlus(x2-x1);
else
    y = x2+SoftPlus(x1-x2);
end
end

function [y] = SoftPlus(x)
if x<-34 % condition for small x
   y = 0;
else
   y = log(1+exp(-x))+x; % numerically stable calculation of log(1+exp(x))
end

end

