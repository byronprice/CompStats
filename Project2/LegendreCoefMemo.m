function [coef] = LegendreCoefMemo(n)
%LegendreCoef.m  Project 2, 1-d WITH MEMOIZATION
%   Code to return coefficients of Legendre's polynomial of order n
%     Implements a computation of the coefficients with effective
%     memoization, though maybe not quite as cleanly as R could do
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
    coef1 = [0,1];
    coef2 = 1;
    
    for ii=2:n
        coef = [0,coef1].*(2*ii-1)-[coef2,0,0].*(ii-1);
        coef = coef./ii;
        
        coef2 = coef1;
        coef1 = coef;
    end
else
    fprintf('Order n must be greater than or equal to 0\n');
    coef = NaN;
end

end