function [X] = GumbelCatSamples(logweights,N)
%GumbelCatSamples.m, Project 4, 2-b
%   draw N samples from a categorical distribution with weights given by
%   exp(logweights) ... uses the Gumbel max trick
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

% Gumbel distribution parameters
sigma = 1;
mu = psi(1); % euler-mascheroni constant

X = zeros(N,1);
for ii=1:N
    % draw K samples from Gumbel distribution  
    G = zeros(K,1);
    for jj=1:K
        Q = rand; % generate uniform random number
        G(jj) = GumbelQuantile(Q,mu,sigma); % generate Gumbel-distributed
                         % random number
    end
    
    [~,ind] = max(G+logweights);
    X(ii) = ind;
end

end
