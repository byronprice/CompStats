function [A,Q,R] = KalmanFilterEM(data)
%KalmanFilterEM.m
%   Fit latent variable (state-space) model for Kalman filter parameters
%    using the EM algorithm. Data is assumed to be driven by an underlying
%    vector autoregressive process
%     e.g. the hidden state space, x, has the following distribution
%          x-(t+1) = A*x-(t)+epsilon
%          where epsilon ~ N(0,Q)
%          
%          this process drives the observed data process,y, by
%          y-(t) = C*x-(t)+nu
%          where nu ~ N(0,R)
%INPUT: data - observed data, input as a matrix, N-by-d, where N is the
%        number of observations and d is the dimensionality of each observation
%OUTPUTS:
%       A - autoregressive(1) process transformation matrix, for the state
%        space x
%       Q - covariance of the state space stochastic process
%       R - covariance of the observation space process
%
%Created: 2018/10/23
% By: Byron Price

% test code with OldFaithfulGeyser dataset
%  GaussMixtureEM(data,2); 
% make sure the data is 272-by-2 columns, with eruption time in minutes
% and wait time to next eruption in minutes

% check size of data
[N,d] = size(data);

dataMean = mean(data,2);

data = data-dataMean;

sigma = cov(data);

% initialize parameters
% generate random covariance matrix for hidden state space noise process
Q = randn(d,d);
Q = Q'*Q;

% estimate of observation noise covariance
R = cov(data);

% could add this as a way to transform x to y, here we assume identity
% C = eye(d);

% generate transformation matrix for vector autoregressive process
A = randn(d,d);
A = A./sum(A(:));

I = eye(d);

Pt1 = eye(d);
Pt = A*Pt1*A'+Q;

K = (Pt+R)\Pt;

PtLag1 = (I-K)*A*Pt;

x = zeros(size(data));

x(2:end,:) = (A*x(1:end-1,:)')';

x(2:end,:) = x(2:end,:)+(K*(data(2:end,:)'-x(2:end,:)'))';

maxIter = 1e3;
tolerance = 1e-6;
for tt=1:maxIter
    % calculate D, E, F
    D = zeros(d,d);E = zeros(d,d);F = zeros(d,d);
    for ii=1:N
        
    end
    
    
    totalPrecision = 0;
    for ii=1:K
       totalPrecision = totalPrecision+sum(abs(muStar{ii}-mu{ii}));
       totalPrecision = totalPrecision+sum(sum(abs(sigmaStar{ii}-sigma{ii})));
    end
    if totalPrecision<=tolerance
        break;
    end

end
end

function [X] = SWEEP(X,k)
%SWEEP.m
%   Code to evaluate the SWEEP operator on row k of matrix X
%Inputs: 
%        X - a matrix, n-by-p
%        k - the row number on which to evaluate the operator, k<=n
%
%Outputs: 
%        X - the matrix after performing the SWEEP operator on row k
%
% Created by: Byron Price
% 2018/10/08

D = X(k,k);
X(k,:) = X(k,:)./D;

[n,~] = size(X);
for ii=[1:k-1,k+1:n]
   B = X(ii,k);
   X(ii,:) = X(ii,:)-B.*X(k,:);
   X(ii,k) = -B/D;
end
X(k,k) = 1/D;
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