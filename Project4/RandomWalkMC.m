function [stopPosition,walkLength,totalPasses] = RandomWalkMC(p)
%RandomWalkMC.m, Project 4, 4-a
%   Monte Carlo sampler for a random walk ... start position is position 1,
%   walker flips a coin with probability of heads p, and if they heads,
%   they move one position to the right (up one number), if tails, they
%   move one position to the left (down one number), stopping conditions
%   are to reach either 0 or 20
%INPUT: p - probability of moving one position to the right (increment)
%OUTPUTS:
%       stopPosition - 0 or 20, the position of the walker at the end of
%          the walk
%       walkLength - whole number greater than 0, number of steps taken
%          throughout the walk
%       totalPasses - the number of times the walker passed position 18
%         (passCondition variable below)
%Created: 2018/10/31
% By: Byron Price

currentPosition = 1;

stopConditions = [0,20];
totalPasses = 0;passCondition = 18; % check to see how many times the walker 
                         % stands on position 18
walkLength = 0;
while currentPosition~=stopConditions(1) && currentPosition~=stopConditions(2)
    stride = 2*binornd(1,p)-1; % step right or left (-1 or +1)
    currentPosition = currentPosition+stride; % take the step
    walkLength = walkLength+1; % add to total steps counter
    if currentPosition==passCondition
        totalPasses = totalPasses+1;
    end
end

stopPosition = currentPosition;
end