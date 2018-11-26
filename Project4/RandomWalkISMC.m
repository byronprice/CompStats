function [stopPosition,walkLength,totalPasses,weight] = RandomWalkISMC(p,pq)
%RandomWalkISMC.m, Project 4, 4-e
%   Monte Carlo sampler for a random walk, with importance sampling
%    ... start position is position 1,
%   walker flips a coin with probability of heads p, and if they heads,
%   they move one position to the right (up one number), if tails, they
%   move one position to the left (down one number), stopping conditions
%   are to reach either 0 or 20
%INPUT: p - probability of moving one position to the right (increment)
%       pq - probability of moving one position to the right for the
%         importance sampler (the actual probability we're using to move)
%OUTPUTS:
%       stopPosition - 0 or 20, the position of the walker at the end of
%          the walk
%       walkLength - whole number greater than 0, number of steps taken
%          throughout the walk
%       totalPasses - the number of times the walker passed position 18
%         (passCondition variable below)
%       weight - importance sampler weight based on the ratio of p to pq
%Created: 2018/10/31
% By: Byron Price

currentPosition = 1;

stopConditions = [0,20];
totalPasses = 0;passCondition = 18; % check to see how many times the walker 
                         % stands on position 18
walkLength = 0;
weight = 1;
while currentPosition~=stopConditions(1) && currentPosition~=stopConditions(2)
    binrnd = binornd(1,pq);
    if binrnd==0
        weight = weight*(1-p)/(1-pq);
    else
        weight = weight*p/pq;
    end
    stride = 2*binrnd-1; % step right or left (-1 or +1)
    currentPosition = currentPosition+stride; % take the step
    walkLength = walkLength+1; % add to total steps counter
    if currentPosition==passCondition
        totalPasses = totalPasses+1;
    end
end

stopPosition = currentPosition;

end