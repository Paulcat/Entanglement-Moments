%Find the minimum distance between ppt and entangled points.
d=[];
for ii=1:length(pptvectors(:,1))
    for mm=1:length(entvectors(:,1))
        d=[d,distance(pptvectors(ii,1), pptvectors(ii,2), pptvectors(ii,3), entvectors(mm,1), entvectors(mm,2), entvectors(mm,3))];
        
    end
end

epsilon=min(d);