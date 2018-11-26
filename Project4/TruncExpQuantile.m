function [X] = TruncExpQuantile(Q,lambda,a,b)
%TruncExpQuantile.m, Project 4, 1-b
%   calculate the inverse CDF of a truncated exponential distribution with 
%   parameters lambda, a, and b
%   for the values provided in Q
%INPUT: Q - values at which to calculate the inverse CDF, from 0 to 1
%       lambda - exponential decay parameter
%       a - left truncation parameter
%       b - right truncation parameter
%OUTPUTS:
%       X - the real-numbered value (greater than 0) corresponding to the
%       quantile Q
%
%Created: 2018/10/30
% By: Byron Price

phi = -1./(exp(-lambda*a)-exp(-lambda*b)); % normalization constant
X = log(Q./phi+exp(-lambda*a))./(-lambda);
end