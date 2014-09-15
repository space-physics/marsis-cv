function [idx distances] = findNearThres(Array,desiredVal)
% this doesn't use the "ismember" technique to then use find, since in this
% project we are using real data that will rarely if ever be an exact match
% with double-precision data
distances = abs(desiredVal-Array);

[~, idx] = min(distances); 

%mDistance = distances(idx);
%nearVal = Array(idx);

end