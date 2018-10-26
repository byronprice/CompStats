function [x,y,Q,R,A] = SimulateKalmanFilterData(N,d)
%SimulateKalmanFilterData.m
%   Simulate data for estimation with a Kalman filter. Observed data is 
%    driven by an underlying vector autoregressive process
%     e.g. the hidden state space, x, has the following distribution
%          x-(t+1) = A*x-(t)+epsilon
%          where epsilon ~ N(0,Q)
%          
%          this process drives the observed data process,y, by
%          y-(t) = C*x-(t)+nu
%          where nu ~ N(0,R)
%INPUT: N - number of observations to simulate
%       d - dimensionality of the observation space
%OUTPUTS:
%       x - the hidden state space data, N-by-d
%       y - observedData, a matrix, N-by-d, of simulated observed data
%       Q - covariance of hidden state noise process
%       R - covariance of observation noise
%       A - transformation matrix for vector autoregressive process
%
%Created: 2018/10/23
% By: Byron Price

y = zeros(N,d);
x = zeros(N,d);

zeroVec = zeros(1,d);

% generate random covariance matrix for hidden state space noise process
Q = [1,0.1,-0.2;0.1,1,0.4;-0.2,0.4,1];

% generate random covariance matrix for observation noise
R = [1,-0.3,0.4;-0.3,1,0.1;0.4,0.1,1];

% generate transformation matrix for vector autoregressive process
A = randn(d,d);
A = A./sum(abs(A(:)));

% A(d,1) = 0;

x(1,:) = mvnrnd(zeroVec,Q);
y(1,:) = x(1,:)+mvnrnd(zeroVec,R);

for ii=2:N
   x(ii,:) = (A*(x(ii-1,:)'))'+mvnrnd(zeroVec,Q);%+sin(ii/N*4*pi);
   y(ii,:) = x(ii,:)+mvnrnd(zeroVec,R);
end

end