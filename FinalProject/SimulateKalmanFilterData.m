function [x,z,Gamma,Sigma,A,C,mu0,V0] = SimulateKalmanFilterData(N,d)
%SimulateKalmanFilterData.m
%   Simulate latent variable (state-space) model for Kalman filter
%    Data is assumed to be driven by an underlying
%    vector autoregressive process
%     e.g. the hidden state space, x, has the following distribution
%          z-(t+1) = A*z-(t)+epsilon
%          where epsilon ~ N(0,Gamma)
%          
%          this process drives the observed data process,y, by
%          x-(t) = C*z-(t)+nu
%          where nu ~ N(0,Sigma)
%INPUT: N - number of observations to simulate
%       d - dimensionality of the observation space
%OUTPUTS:
%       x - "observed" data, d-by-N
%       z - the hidden states, d-by-N
%       A - auto-regressive(1) process transformation matrix, for the state
%        space z
%       C - transformation from hidden state space to observations
%       Gamma - covariance of the state space stochastic process
%       Sigma - covariance of the observation space process
%       mu0 - state space initial conditions mean
%       V0 - state-space covariance initial condition
%
%Created: 2018/12/11
% By: Byron Price

z = zeros(d,N);
x = zeros(d,N);

zeroVec = zeros(1,d);

% generate random covariance matrix for hidden state space noise process
Gamma = [1,0.1,-0.2;0.1,1,0.4;-0.2,0.4,1]./5;

% generate random covariance matrix for observation noise
Sigma = [1,-0.3,0.4;-0.3,1,0.1;0.4,0.1,1];

% generate transformation matrix for vector autoregressive process
% A = randn(d,d);
% A = A./sum(abs(A(:)));

A = [1,-0.3,0.4;-0.3,1,0.1;0.4,0.1,1];
A = A./sum(abs(A));

C = eye(d);
% A(d,1) = 0;
mu0 = zeroVec';V0 = Gamma;

z(:,1) = mvnrnd(mu0',V0)';
x(:,1) = C*z(:,1)+mvnrnd(zeroVec,Sigma)';

for ii=2:N
   z(:,ii) = A*z(:,ii-1)+mvnrnd(zeroVec,Gamma)';%+sin(ii/N*4*pi);
   x(:,ii) = C*z(:,ii)+mvnrnd(zeroVec,Sigma)';
end

end