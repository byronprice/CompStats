function [X] = CategoricalSamples(logweights,N)
%CategoricalSamples.m, Project 4, 2-a
%   draw N samples from a categorical distribution with weights given by
%   exp(logweights) ...
%INPUT: logweights - the log of the weights of the distribution, the length
%        of the vector is the number K of categories
%       N - number of random samples to draw
%OUTPUTS:
%       X - the N random samples, with values from 1 to K indicating which
%        category the sample belongs to
%
%Created: 2018/10/30
% By: Byron Price
K = length(logweights);
logweights = logweights(:);
logpi = logweights-LogSum(logweights,K); % log probabilities for each category

logci = zeros(K,1);

for ii=1:K
    logci(ii) = LogSum(logpi(1:ii),ii); % log CDF
end

X = zeros(N,1);
for ii=1:N
    logU = log(rand);
    check = min(logci-logU,0);
    ind = find(check==0,1,'first');
    X(ii) = ind;
end

end

function [summation] = LogSum(vector,vectorLen)
% calculate log of sum of exponentials
%   ie log(sum[ exp(vector) ] )
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