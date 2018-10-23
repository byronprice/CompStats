function [Q,nu,alpha] = Q_Vandermonde2(x,d)
%Q_Vandermonde2.m  Project 1, 3-c
%   Calculate the Q orthogonal factor for Vandermonde matrix of vector x 
%    with degree d ... calculated using a recursive relation
%    each column of the Vandermonde matrix is a term in the geometric
%    progression of a (first column is x^0, second is x^1, etc.)
%Inputs:
%        x - an n-by-1 traditional (column) vector
%        d - degree of the final matrix
%Outputs:
%        Q - Q orthogonal factor for the Vandermonde matrix of vector x
%        with degree d
%        nu - d+2-by-1 vector of variables used in recursion
%        alpha - d+1-by-1 vector of variables used in recursion
%
% Created by: Byron Price
% 2018/09/21
x = x(:);n = length(x); % make x a column vector, get its length

diagX = diag(x); % calculate sufficient "stat", diagonal matrix with x on diag

U = ones(n,d+1); % initialize U

nu = zeros(d+2,1);nu(1) = 1; % initialize nu
alpha = zeros(d+1,1); % initialize alpha

% get second column of U
nu(2) = GetNu(U(:,1)); 
alpha(1) = GetAlpha(U(:,1),diagX);

U(:,2) = x-alpha(1).*ones(n,1);

for ii=2:d % iterate
   nu(ii+1) = GetNu(U(:,ii));
   alpha(ii) = GetAlpha(U(:,ii),diagX);
   for jj=1:n
       U(jj,ii+1) = (x(jj)-alpha(ii))*U(jj,ii)-(nu(ii+1)/nu(ii))*U(jj,ii-1);
   end
end
alpha(ii+1) = GetAlpha(U(:,ii+1),diagX);
nu(ii+2) = GetNu(U(:,ii+1));

Q = normalizeCols(U); % column normalize U
end

function [alpha] = GetAlpha(u_i,diagX)
alpha = (u_i'*diagX*u_i)/(u_i'*u_i);
end

function [nu] = GetNu(u_i)
nu = u_i'*u_i;
end


