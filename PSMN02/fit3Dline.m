%% fit3Dline
function [xyz0,direction] = fit3Dline(XYZ)

if max(max(max(isnan(XYZ)))) ==0
    [xyz0,direction] = fit3Dline_nonan(XYZ);
else
    [P V] = arrayfun(@(I)(fit3Dline_nan(XYZ(:,:,I))),1:size(XYZ,3),'UniformOutput',false);
    xyz0 = (cell2mat(P'));
    direction = (cell2mat(V'));
    
    xyz0(isnan(xyz0)) = [];
    direction(isnan(direction)) = [];
end

end

