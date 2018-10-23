function [eps] = epsilon(start)
%epsilon.m  Project 1, 2-a
%   Calculate the machine epsilon for the computer that is running this
%   code
%
%Inputs: 
%        eps - a starting number (defaults to 1)
%Outputs: 
%        eps - the machine epsilon
%
% Created by: Byron Price
% 2018/09/21

if nargin<1
    eps = 1;
else
    eps = start;
end

while (1+eps/2)~=1 % machine precision definition reached when
                 % 1+eps/2 == 1
    eps = eps/2;
end

end

