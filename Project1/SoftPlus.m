function [y] = SoftPlus(x)
%SoftPlus.m  Project 1, 2-b
%   Calculate the softplus function of the input x
%    softplus(x) = log(1+exp(x))
%    it is a smooth approximation to the half-wave rectifier function
%     f(x) = max(0,x)
%
%Inputs: 
%        x - any real number 
%Ouputs:
%        y - the output of the softplus function
%
% Created by: Byron Price
% 2018/09/21

if x<-34 % condition for small x
   y = 0;
else
   y = log(1+exp(-x))+x; % numerically stable calculation of log(1+exp(x))
end

end

