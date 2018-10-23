function [hornerFun] = myHorner(coef)
%myHorner.m  Project 2, 1-a
%   Code to create function that will evaluate polynomial with Horner's scheme
%    p(x) = a0+a1*x+a2*x^2+...+an*x^n
%Inputs: 
%        coef - a vector of polynomial coefficients [a0,a1,...,an]
%
%Outputs: 
%        hornerFun - function handle that can be used to evaluate the 
%         polynomial defined by coef at a point x
%         e.g. hornerFun(x) = p(x)
%
% Created by: Byron Price
% 2018/10/08

hornerFun = @(x) funcDef(x,coef);

end

function [y] = funcDef(x,coef)

y = zeros(size(x));

for ii=length(coef):-1:1
    y = y.*x+coef(ii);
end
end
