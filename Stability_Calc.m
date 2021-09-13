function [stability] = calculateStability(sData, skipTime, skipTime2)

findHet = sData(:,2) == 1;
sData = sData(~findHet,[1 3]);
findTimes = sData(:,2) >= skipTime & sData(:,2) <= skipTime2;
sData = sData(findTimes,:) ;
[x,y] = sort(sData(:,1))   ;
sData = sData(y,:);

findMaxTime = sData(:,2) >= min(sData(:,2)) + (max(sData(:,2))-min(sData(:,2)))/2;

if(isempty(sData))
    stability = 0;
    return;
end

AMD_1 = AMDv4(sData(~findMaxTime,:));
AMD_2 = AMDv4(sData(findMaxTime,:));

stability = AMDv4_Similarity_c(AMD_1, AMD_2, min([min(AMD_1(1,:)) min(AMD_2(1,:))]), max([max(AMD_1(1,:)) max(AMD_2(1,:))]));
stability = stability(2);

end


