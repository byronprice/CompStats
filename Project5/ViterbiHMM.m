function [logProbPath,states] = ViterbiHMM(P,EmissionDist,start,emission)
%ForwardHMM.m   Project 5, 1-b
%   Implements the Viterbi algorithm
%    given a Hidden Markov model with state transition probabilities
%   given by P and assuming the emission distribution given the state is a
%   normal random variable with a unique mean and variance for each state,
%   find the most likely sequence of hidden states
%INPUTS:
%        P - transition probability matrix [number of states by number of
%           states]
%        Emissions - emission matrix (mean and variance of normal
%           distribution for each state), [number of states by 2]
%        start - state that chain starts in
%        emission - observed data
%OUTPUTS:
%        logProbPath - probability of maximum-likelihood path, given the
%           observed data
%        states - the most likely sequence of hidden states
% Byron Price, 2018/11/25
N = size(emission,1);

logP = log(P);
numStates = size(P,1);

logNormPDF = @(x,mu,sigmasquare) -0.5*log(2*pi*sigmasquare)-((x-mu).^2)./(2*sigmasquare);

V = zeros(N,numStates);
B = zeros(N,numStates);

for jj=1:numStates
    logygivenx = logNormPDF(emission(1),EmissionDist(jj,1),EmissionDist(jj,2));
    logxgivenx = logP(start,jj);
    
    V(1,jj) = logygivenx+logxgivenx;
    B(1,jj) = 0;
end
for ii=2:N
    for jj=1:numStates
        logygivenx = logNormPDF(emission(ii),EmissionDist(jj,1),EmissionDist(jj,2));
        
        logVec = zeros(numStates,1);
%        logVec2 = zeros(numStates,1);
        for kk=1:numStates
            logxgivenx = logP(kk,jj);
            logVec(kk) = logxgivenx+logygivenx+V(ii-1,kk);
%             logVec2(kk) = logxgivenx+V(ii-1,kk);
        end
        [val,ind] = max(logVec);
        V(ii,jj) = val;
%         [~,ind] = max(logVec2);
        B(ii,jj) = ind;
    end
end
[val,ind] = max(V(end,:));
logProbPath = val;

% backtrace
states = zeros(N,1);

states(N) = ind;
for ii=N-1:-1:1
    states(ii) = B(ii+1,states(ii+1));
end


end