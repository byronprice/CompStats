function [coef] = LegendreCoef(n)
%LegendreCoef.m  Project 2, 1-d NO MEMOIZATION
%   Code to return coefficients of Legendre's polynomial of order n
%     Implements a recursive computation of the coefficients without 
%      memoization
%Inputs: 
%        n - the order of the polynomial who's coefficients you wish to return
%
%Outputs: 
%        coef - a vector of the Legendre polynomial coefficients
%          if p(x) = a0+a1*x+a2*x^2+...+an*x^n , then 
%          coef = [a0,a1,a2,...,an]
%
% Created by: Byron Price
% 2018/10/08

if n==0
    coef = 1;
elseif n==1
    coef = [0,1];
elseif n>1
    coef1 = LegendreCoef(n-1);
    coef2 = LegendreCoef(n-2);
    coef = [0,coef1].*(2*n-1)-[coef2,0,0].*(n-1);
    coef = coef./n;
else
    fprintf('Order n must be greater than or equal to 0\n');
    coef = NaN;
end

end