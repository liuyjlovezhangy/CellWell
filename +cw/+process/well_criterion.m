function [good_well,o] = well_criterion(object)
    o.min_well_area = 2000;
    o.max_well_area = 5000;
    
    if ~exist('object','var')
        good_well = [];
        
        return
    end
    
    if object.Area <o. min_well_area || object.Area > o.max_well_area
        good_well = 0;
    else
        good_well = 1;
    end
    
end