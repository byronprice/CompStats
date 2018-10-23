function [Q] = Q_Vandermonde(x,d)
%Q_Vandermonde.m  Project 1, 3-b
%   Calculate the Q orthogonal factor for Vandermonde matrix of vector a 
%    with degree d, calculated based on Gram-Schmidt orthogonalization
%    each column of the Vandermonde matrix is a term in the geometric
%    progression of a (first column is a^0, second is a^1, etc.)
%Inputs:
%        a - an n-by-1 traditional (column) vector
%        d - degree of the final matrix
%Outputs:
%        Q - Q orthogonal factor for the Vandermonde matrix of vector a
%        with degree d
%
% Created by: Byron Price
% 2018/09/21
x = x(:);n = length(x); % ensure x is column vector, get it's length

U = ones(n,d+1); % initialize U

for ii=2:d+1 % iterate
   vandermondeVec = x.^(ii-1); % create vandermonde vector, a_i
   U(:,ii) = vandermondeVec; % initialize u_i with a_i
   for jj=1:ii-1
      U(:,ii) = U(:,ii)-proj(vandermondeVec,U(:,jj)); % sum over projections 
                          % to create orthogonal set
   end
end

Q = normalizeCols(U); % column normalize
end

