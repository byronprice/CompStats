function [root] = NewtonRaphson(evalFunc,evalFuncDeriv,startVal)
%NewtonRaphson.m  Project 2, 1-c
%   Code to evaluate roots of a function using the Newton-Raphson method
%Inputs: 
%        evalFunc - function handle for the function who's roots you wish
%         to evaluate
%        evalFuncDeriv - function handle to the derivative of the function
%         evalFunc
%
%Outputs: 
%        root - the root (zero) of the function closest to the point given
%         by startVal
%
% Created by: Byron Price
% 2018/10/08

maxIter = 1e4;
tolerance = 1e-6;

root = startVal;
for ii=1:maxIter
   rootprime = root-evalFunc(root)/evalFuncDeriv(root);
   if abs(evalFunc(rootprime)-evalFunc(root))<tolerance
       root = rootprime;
       break;
   end
   root = rootprime;
end

fprintf('Number of iterations: %d\n',ii);
end
