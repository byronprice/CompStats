function [logProbData,logAlpha] = ForwardHMM(P,EmissionDist,start,emission)
%ForwardHMM.m   Project 5, 1-a
%   Implements the forward algorithm
%    given a Hidden Markov model with state transition probabilities
%   given by P and assuming the emission distribution given the state is a
%   normal random variable with a unique mean and variance for each state,
%   find the probability of the data P(y)
%INPUTS:
%        P - transition probability matrix [number of states by number of
%           states]
%        Emissions - emission matrix (mean and variance of normal
%           distribution for each state), [number of states by 2]
%        start - state that chain starts in
%        emission - observed data
%OUTPUTS:
%        logProbData - probability of the observed data, under the model
%        logAlpha - log of the probability of being in a given state at a
%          given moment
% Byron Price, 2018/11/25
N = size(emission,1);

logP = log(P);
numStates = size(P,1);
prevAlpha = -Inf*ones(1,numStates);prevAlpha(start) = 0;

logAlpha = zeros(N,numStates);
logNormPDF = @(x,mu,sigmasquare) -0.5*log(2*pi*sigmasquare)-((x-mu).^2)./(2*sigmasquare);
for ii=1:N
    for jj=1:numStates
        logygivenx = logNormPDF(emission(ii),EmissionDist(jj,1),EmissionDist(jj,2));
        
        logVec = zeros(numStates,1);
        for kk=1:numStates
            logxgivenx = logP(kk,jj);
            logVec(kk) = logxgivenx+prevAlpha(kk);
        end
        logAlpha(ii,jj) = logygivenx+LogSum(logVec,numStates);
%         vec = zeros(numStates,1);
%         for kk=1:numStates
%             xgivenx = P(kk,jj);
%             vec(kk) = xgivenx*exp(prevAlpha(kk));
%         end
%         logAlpha(ii,jj) = logygivenx+log(sum(vec));
    end
    prevAlpha = logAlpha(ii,:);
end
logProbData = LogSum(logAlpha(N,:),numStates);

end