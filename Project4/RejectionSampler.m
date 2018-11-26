function [X] = RejectionSampler(N)
%RejectionSampler.m, Project 4, 3-c
%   Rejection sampler to generate random samples from the standard normal
%   distribution
%INPUT: N - number of samples to generate
%OUTPUTS:
%       X - the random samples, from the standard normal
%
%Created: 2018/10/30
% By: Byron Price

X = zeros(N,1);

% sample categorical for mixture of truncated exp's
pii = [0.3799,0.4803,0.1398]; % mixture weights for each component
logweights = log(pii);

lambda1 = 1;a1 = 1/2;b1 = 3/2; % truncated exp 1
lambda2 = 2;a2 = 3/2;b2 = Inf; % truncated exp 2

M = 0.525; % M parameter for rejection sampling
% reject = [];
for ii=1:N
    U = 2;
    check = 1;
    while U>=check
        mixture = GumbelCatSamples(logweights,1);
        if mixture==1
            xstar = rand*0.5;
            q = pii(1)*2;
        elseif mixture==2
            Q = rand;
            xstar = TruncExpQuantile(Q,lambda1,a1,b1);
            q = pii(2)*(exp(3/2)/(exp(1)-1))*exp(-lambda1*xstar);
        elseif mixture==3
            Q = rand;
            xstar = TruncExpQuantile(Q,lambda2,a2,b2);
            q = pii(3)*2*exp(1)*exp(-lambda2*xstar);
        end
        U = rand;
        p = normpdf(xstar,0,1);
        check = p/(M*q);
%         reject = [reject;U<check];
    end
    sign = binornd(1,0.5)*2-1; % get either positive or negative values
    X(ii) = xstar*sign;
end
% reject = mean(reject);
end