function [states,steadyStateProb] = SecondOrderMarkovChain(T,N,initProb)
% SecondOrderMarkovChain.m
%  T - transition probability matrix
%  N - length of the chain
%  (optional) initProb - initial probabilities (defaults to uniform)

% example transition matrix
% T = [0,0,0,0,0;1/4,0,1/4,1/4,1/4;0.1/3,0.1/3,0,0.1/3,0.9;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;0,0,0,0,0;0.9,0.1/3,0,0.1/3,0.1/3;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0,0,0,0,0;1/4,1/4,1/4,0,1/4;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0.1/3,0.9,0,0.1/3,0.1/3;0,0,0,0,0;1/4,1/4,1/4,1/4,0];
% T = [T;0,1/4,1/4,1/4,1/4;1/4,0,1/4,1/4,1/4;0.1/3,0.1/3,0,0.9,0.1/3;1/4,1/4,1/4,0,1/4;0,0,0,0,0];

numStates = size(T,2);

states = zeros(N,1);

if nargin<3
    initProb = ones(1,numStates)./numStates;
end

states(1) = find(mnrnd(1,initProb));
states(2) = find(mnrnd(1,initProb));

while states(1)==states(2)
    states(2) = find(mnrnd(1,initProb));
end

stateMat = zeros(numStates*numStates,2);
index = 1;
for ii=1:numStates
    for jj=1:numStates
        stateMat(index,1) = ii;stateMat(index,2) = jj;
        index = index+1;
    end
end

for ii=3:N
   oneBack = states(ii-1);
   twoBack = states(ii-2);
   index = intersect(find(stateMat(:,2)==oneBack),find(stateMat(:,1)==twoBack));
   transitionProb = T(index,:);
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