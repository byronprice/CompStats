function [states,steadyStateProb] = SecondOrderMarkovChainTwo(T2,T1,N,initProb)
% SecondOrderMarkovChain.m
%  T2 - two back transition probability matrix, T1 - one back transition probability
%   matrix
%  N - length of the chain
%  (optional) initProb - initial probabilities (defaults to uniform)

% example transition matrix
% T = [0,0,0,0,0;1/4,0,1/4,1/4,1/4;0.1/3,0.1/3,0,0.1/3,0.9;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;0,0,0,0,0;0.9,0.1/3,0,0.1/3,0.1/3;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0,0,0,0,0;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0.1/3,0.9,0,0.1/3,0.1/3;0,0,0,0,0;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0.1/3,0.1/3,0,0.9,0.1/3;1/4,1/4,1/4,0,1/4;0,0,0,0,0];

numStates = size(T1,2);

states = zeros(N,1);

if nargin<4
    initProb = ones(1,numStates)./numStates;
end

states(1) = find(mnrnd(1,initProb));
states(2) = find(mnrnd(1,initProb));

while states(1)==states(2)
    states(2) = find(mnrnd(1,initProb));
end

for ii=3:N
   oneBack = states(ii-1);
   twoBack = states(ii-2);
   
   transitionProb = T1(oneBack,:).*T2(twoBack,:);
   transitionProb = transitionProb./sum(transitionProb);
   states(ii) = find(mnrnd(1,transitionProb));
end

steadyStateProb = ones(1,numStates)./numStates;
% [~,D,W] = eig(T);tolerance = 1e-6;
% for ii=1:numStates
%     if abs(D(ii,ii)-1)<tolerance
%        steadyStateProb = W(:,ii)'./sum(W(:,ii));
%     end
% end

end