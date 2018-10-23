function [vanMatrix] = vandermonde(a,d)
%vandermonde.m  Project 1, 1-e
%   Calculate the Vandermonde matrix of vector a with degree d
%    each column of the Vandermonde matrix is a term in the geometric
%    progression of a (first column is a^0, second is a^1, etc.)
%Inputs:
%        a - an n-by-1 traditional (column) vector
%        d - degree of the final matrix
%Outputs:
%        vanMatrix - vandermonde matrix of a
%
% Created by: Byron Price
% 2018/09/21
a = a(:); % ensure a is a column vector
n = length(a);

vanMatrix = ones(n,d+1); % initialize matrix

for ii=2:d+1
   vanMatrix(:,ii) = a.^(ii-1); % take geometric progression
end
end

