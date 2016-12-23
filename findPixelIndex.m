function RowColIndex = findPixelIndex(AltInterp,freqLinInterp,pickedPoint)

% if size(pickedPoint,2)>2, N = 1;    %this was from a mouse click
% else N = size(pickedPoint,1);       %this was from preset coordinates
% end
% 
% RowColIndex = nan(N,2); %preinitialize
% 
% for i = 1:N
RowColIndex(:,1) = round(interp1(AltInterp,1:length(AltInterp),pickedPoint(:,2)));

RowColIndex(:,2) = round(interp1(freqLinInterp,1:length(freqLinInterp),pickedPoint(:,1)));
%end

end