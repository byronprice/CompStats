function [X] = GumbelQuantile(Q,mu,sigma)
%GumbelQuantile.m, Project 4, 1-a
%   calculate the inverse CDF of a Gumbel distribution with parameters 
%   mu and sigma for the values provided in Q
%INPUT: Q - values at which to calculate the Gumbel inverse CDF, a number
%         from 0 to 1
%       mu - mean parameter (location)
%       sigma - variance parameter (scale)
%OUTPUTS:
%       X - the real-numbered value corresponding to a given quantile
%
%Created: 2018/10/30
% By: Byron Price

X = mu-sigma.*log(-log(Q));


end