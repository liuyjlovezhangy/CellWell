function linked_object_cell = assign_tracks_to_objects( tracks, object_ids, objects_cell )
    % linked_object_cell: structure 1 x num unique cells in well over time
    %   start_frame
    %   end_frame
    %   cell_props (1 x end_frame - start_frame + 1)
    %   track

    warning('assign_tracks_to_objects: linked_object_cell addition needs to be vectorized.')
    
    num_tracks = size(object_ids,1);
    num_cells = num_tracks;
    
    for cell_idx = 1:num_cells
        cur_object_ids = object_ids(cell_idx,:);
        start_frame = find(~isnan(cur_object_ids),1,'first');
        end_frame = find(~isnan(cur_object_ids),1,'last');
        cur_track = tracks{cell_idx};
        
        linked_object_cell(cell_idx).start_frame = start_frame;
        linked_object_cell(cell_idx).end_frame = end_frame;
        linked_object_cell(cell_idx).track = cur_track;
        
        %%% Add all object properties entries to this cell structure entry
        
        cell_props = [];
        
        for frame_idx = start_frame:end_frame
            cell_props = [cell_props, objects_cell{frame_idx}(cur_object_ids(frame_idx))];
        end
        
        linked_object_cell(cell_idx).cell_props = cell_props;
        
    end
    
end

