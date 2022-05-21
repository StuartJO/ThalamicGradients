function [xrange,yrange,zrange] = getLimitsRange(startRange,endRange,timepoints,startStick,endStick)

if nargin < 5
    startStick = [];
end

if nargin < 5
    endStick = [];
end

if ~isempty(startStick)
for i = 1:3
range_start{i} = repmat(startRange(i,:),startStick,1);
end
timepoints = timepoints-startStick;
else
    for i = 1:3
   range_start{i} = [];
    end
end
if ~isempty(endStick)
for i = 1:3
range_end{i} = repmat(endRange(i,:),endStick,1);
end
timepoints = timepoints-endStick;
else
    for i = 1:3
   range_end{i} = [];
    end
end

for i = 1:3
range{i} = [range_start{i}; linspace(startRange(i,1),endRange(i,1),timepoints)' linspace(startRange(i,2),endRange(i,2),timepoints)';range_end{i}];
end

xrange = range{1};
yrange = range{2};
zrange = range{3};