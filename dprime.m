function [dPrime, c] = dprime(HIT,FA,N)
%function [dPrime c] = dprime(HIT,FA,N)
% This function calculates the criterion c and information d' values.
% use as: dprime(HIT,FA,N)
% INPUTS
% HIT - Insert the number (not proportion) of correct trials
% FA - Insert the number of false alarm trials
% N - number of trials for Hit and FA - for calculating proportions
% OUTPUT
% dprime - dPrime value
% c      - bias (if 0 - no bias, >0 - biased toward 'yes', <0 - biased toward 'no'


% adjustments if all correct/wrong
if HIT == N || FA == 0
   HIT = N-1;
   FA = 1/N;
elseif FA == N || HIT == 0
   FA = N-1;
   HIT = 1/N;    
end

zHit = norminv(HIT/N);
zFA = norminv(FA/N);

% calculate d-prime
dPrime = zHit - zFA;

% calculate bias / criterion
c = (zHit + zFA)/2;
end
