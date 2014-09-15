function B = BackgroundModes(wgt,sigma,BackgroundThres)
% by Michael Hirsch
% March 11, 2012
% based off of Stauffer and Grimson GMM in 1999 paper.
%% background detection
% How many distribution are necessary to represent the background? A simple
% unimodal model works for a completely constant pixel. More complicated
% backgrounds (ripples in water, waving flag, etc.) require additional
% distribution (2, 3, 4, etc. modal Gaussian distribution)

%functionally, here are the steps
% (observe, we sort by weights/sigma, but cumsum by weights alone)
% 1) sort J = weights/sigma 
% 2) starting with the biggest weight, keep adding the weights (not weights/sigma)
% in order from greatest->least, until BackgroundThres is exceeded-- this is "B",
% the number of modes in the background distributions

B = 0; wSum = 0;
[~, iJ] = sort(wgt./sigma,2,'descend');
while wSum <= BackgroundThres
B = B+1;
wSum = wSum + wgt(iJ(B)); %wgt(iJ(B)) picks the next weight to add, in descending order of wgt/sigma
end
end