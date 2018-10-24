function [x,P,K] = RealTimeKalmanFilter(y,Q,R,A)
%RealTimeKalmanFilterData.m
%   Use Kalman filter to predict unobserved state space. Observed data is 
%    driven by an underlying vector autoregressive process
%     e.g. the hidden state space, x, has the following distribution
%          x-(t+1) = A*x-(t)+epsilon
%          where epsilon ~ N(0,Q)
%          
%          this process drives the observed data process,y, by
%          y-(t) = C*x-(t)+nu
%          where nu ~ N(0,R)
%INPUT: y - observed measurement data, N-by-d, observations by dimensions
%       Q - state space noise covariance
%       R - observation space noise covariance
%       A - vector autoregressive process transformation matrix
%OUTPUTS:
%       x - the hidden state space data, N-by-d
%       P - expected error covariance matrix for prediction of x
%       K - Kalman gain
%
%Created: 2018/10/23
% By: Byron Price

[N,d] = size(y);

meanY = mean(y);

y = y-meanY;

y = y';

x = zeros(d,N);

I = eye(d);

P = eye(d);
Phat = A*P*A'+Q;
K = (Phat+R)\Phat;
P = (I-K)*A*Phat;

%tolerance = 0;
for ii=2:N
    xmeas = A*x(:,ii);
    Phat = A*P*A'+Q;
    newK = (Phat+R)\Phat;
    P = (I-newK)*A*Phat;
    
    x(:,ii) = xmeas+newK*(y(:,ii)-xmeas);

%     if sum(sum(abs(K-newK)))==tolerance
%         break;
%     end
    
    K = newK;
end

x(:,1) = K*y(:,1);

for ii=2:N
    xmeas = A*x(:,ii);
    x(:,ii) = xmeas+K*(y(:,ii)-xmeas);
end

x = x'+meanY;
end