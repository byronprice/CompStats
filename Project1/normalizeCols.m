function [Matrix] = normalizeCols(Matrix)
%normalizeCols.m  Project 1, 1-c
%   Normalize each column of the input Matrix by its column L2 norm
%Inputs: 
%        Matrix - an n-by-m matrix
%Outputs:
%        Matrix - the input matrix with each column normalized by its L2 norm
%
% Created by: Byron Price
% 2018/09/21
[~,m] = size(Matrix);

% iterate through each column
for ii=1:m
   L2norm = norm2(Matrix(:,ii)); % get L2 norm
   Matrix(:,ii) = Matrix(:,ii)./L2norm; % divde column by its L2 norm
end

end

