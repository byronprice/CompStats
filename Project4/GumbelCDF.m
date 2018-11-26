function [F] = GumbelCDF(x,mu,sigma)
%GumbelCDF.m, Project 4, 1-a
%   calculate the CDF of a Gumbel distribution with parameters mu and sigma
%   for the values provided in x
%INPUT: x - values at which to calculate the Gumbel CDF, any real number
%       mu - mean parameter (location)
%       sigma - variance parameter (scale)
%OUTPUTS:
%       F - the CDF (e.g. Probability(X<=x))
%
%Created: 2018/10/30
% By: Byron Price

F = exp(-exp(-(x-mu)./sigma));


end