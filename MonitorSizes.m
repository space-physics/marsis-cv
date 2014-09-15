function MS = MonitorSizes()
MS = get(0,'MonitorPositions');
%choose to use the largest monitor for displays
if size(MS,1)>1
if MS(2,3)>MS(1,3) %then swap rows
    MS = flipud(MS);
end
end
end