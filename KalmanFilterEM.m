function [A,Q,R,x] = KalmanFilterEM(data)
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
%       x - estimates of the hidden states, N-by-d
%
%Created: 2018/10/23
% By: Byron Price

% test code with OldFaithfulGeyser dataset
%  GaussMixtureEM(data,2); 
% make sure the data is 272-by-2 columns, with eruption time in minutes
% and wait time to next eruption in minutes

% check size of data
[N,d] = size(data);

dataMean = mean(data,1);

data = data-dataMean;

data = data';

% sigma = cov(data');

% initialize parameters
% generate random covariance matrix for hidden state space noise process
Q = [1,0.1,-0.2;0.1,1,0.4;-0.2,0.4,1];

% estimate of observation noise covariance
R = [1,-0.3,0.4;-0.3,1,0.1;0.4,0.1,1];

% could add this as a way to transform x to y, here we assume identity
% C = eye(d);

% generate transformation matrix for vector autoregressive process
tmp1 = zeros(d,d);
tmp2 = zeros(d,d);
for ii=2:N
    tmp1 = tmp1+data(:,ii-1)*data(:,ii)';
    tmp2 = tmp2+data(:,ii-1)*data(:,ii-1)';
end
A = (tmp2\tmp1)';

x = zeros(d,N+1);

maxIter = 1e3;
tolerance = 1e-4;
for tt=1:maxIter
    % E step, get expected hidden state estimates
    Rinv = R;
    Qinv = Q;
    for ii=1:d
        Rinv = SWEEP(Rinv,ii); % inverse of R and Q now saved
        Qinv = SWEEP(Qinv,ii);
    end
%     x(:,1) = mvnrnd(zeros(d,1),Q)';
    for ii=2:N+1
        x(:,ii) = (Rinv+Qinv)\(Rinv*data(:,ii-1)+Qinv*A*x(:,ii-1));
    end
    
    % M step, maximize expected log-likelihood
    tmp1 = zeros(size(A));
    tmp2 = zeros(size(A));
    for ii=2:N+1
        tmp1 = tmp1+x(:,ii-1)*x(:,ii)';
        tmp2 = tmp2+x(:,ii-1)*x(:,ii-1)';
    end
    Ahat = (tmp2\tmp1)';
    Qhat = zeros(size(Q));
    Rhat = zeros(size(R));
    
    for ii=2:N+1
        Qhat = Qhat+(x(:,ii)-Ahat*x(:,ii-1))*(x(:,ii)-Ahat*x(:,ii-1))';
        Rhat = Rhat+(data(:,ii-1)-x(:,ii))*(data(:,ii-1)-x(:,ii))';
    end
    Qhat = Qhat./N;
    Rhat = Rhat./N;
    
    
    totalPrecision = abs(A-Ahat)+abs(Q-Qhat)+abs(R-Rhat);
    if sum(totalPrecision(:))<=tolerance
        break;
    end
    A = Ahat;
    Q = Qhat;
    R = Rhat;
end

x = x(:,2:end)'+dataMean;
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