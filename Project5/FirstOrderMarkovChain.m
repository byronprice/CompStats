function [states,steadyStateProb] = FirstOrderMarkovChain(T,N,initProb)
% FirstOrderMarkovChain.m
%  T - transition probability matrix
%  N - length of the chain
%  (optional) initProb - initial probabilities (defaults to uniform)

% example transition matrices
% T = [0,0.8,0.1,0.1;0.1,0,0.8,0.1;0.1,0.1,0,0.8;0.8,0.1,0.1,0];
% val3 = (1-0.25-0.1/3)/2;
% T = [0,0.9,0.1/3,0.1/3,0.1/3;0.25,0,0.25,0.25,0.25;0.25,0.1/3,0,val3,val3;0.25,0.1/3,val3,0,val3;0.25,0.1/3,val3,val3,0];

numStates = size(T,1);

states = zeros(N,1);

if nargin<3
    initProb = ones(1,numStates)./numStates;
end

states(1) = find(mnrnd(1,initProb));

for ii=2:N
   oneBack = states(ii-1);
   transitionProb = T(oneBack,:);
   states(ii) = find(mnrnd(1,transitionProb));
end

[~,D,W] = eig(T);tolerance = 1e-6;
for ii=1:numStates
    if abs(D(ii,ii)-1)<tolerance
       steadyStateProb = W(:,ii)'./sum(W(:,ii));
    end
end

end