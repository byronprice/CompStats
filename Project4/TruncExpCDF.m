function [F,X] = TruncExpCDF(x,lambda,a,b)
%TruncExpCDF.m, Project 4, 1-b
%   calculate the CDF of a truncated exponential distribution with 
%   parameters lambda, a, and b
%   for the values provided in x
%INPUT: x - values at which to calculate the CDF, any real number >=0
%       lambda - exponential decay parameter
%       a - left truncation parameter
%       b - right truncation parameter
%OUTPUTS:
%       F - the CDF (e.g. Probability(X<=x))
%
%Created: 2018/10/30
% By: Byron Price
x = x(x<b);

phi = -1./(exp(-lambda*a)-exp(-lambda*b)); % normalization constant

F = (exp(-lambda.*x)-exp(-lambda*a)).*phi;

I = x>=a;
F = F.*I;
X = x;
end