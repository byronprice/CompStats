% cars.m
%  use cars dataset to calculate ML estimate for stopping distance as a
%  function of speed
%   going to make a Vandermonde design matrix, such that the stopping
%   distance will be polynomial function of the speed

M = dlmread('cars.csv',',',1,1);

x = M(:,1);y = M(:,2);

d = 3;

X = vandermonde(x,d); % get Vandermonde design matrix

% compute gamma with Q from the QR decomposition
Q = Q_Vandermonde(x,d);

p = size(Q,2);

gamma = Q'*y;

[gamma2,~,stats] = glmfit(Q,y,'normal','constant','off');

stdev = std(y-Q*gamma);
confInt = zeros(p,2);
keepParams = ones(p,1);
for ii=1:p
    confInt(ii,1) = gamma(ii)-1.96*stdev;confInt(ii,2) = gamma(ii)+1.96*stdev;
    if confInt(ii,1)<0 && confInt(ii,2)>0
        keepParams(ii) = 0;
        gamma(ii) = 0;
        confInt(ii,:) = 0;
    end
end

% get inverse of R, so that we can find beta
Rinv = X\Q;

beta_hat = Rinv*gamma;

beta_confInt = Rinv*confInt;

% NOW TRY TO GET BETA THE SECOND WAY

[beta,~,stats2] = glmfit(X,y,'normal','constant','off');

[beta2,corrMat] = MLestimate(X,y);