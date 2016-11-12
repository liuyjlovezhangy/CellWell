function linked_object_cell = assign_tracks_to_objects( tracks, object_ids, objects_cell )
    % linked_object_cell: structure 1 x num unique cells in well over time
    %   start_frame
    %   end_frame
    %   regionprops (1 x end_frame - start_frame + 1)
    %   track

    num_tracks = size(object_ids,1);
    num_frames = size(object_ids,2);
    num_cells = num_tracks;
    
    object_ids
    
    linked_object_cell = [];
    
    for cell_idx = 1:num_cells
        cur_object_ids = object_idcs(cell_idx,:);
        
        linked_object_cell(cell_idx).start_frame = find(~isnan(cur_object_ids),1,'first');
        linked_object_cell(cell_idx).end_frame = find(~isnan(cur_object_ids),1,'last');
        linked_object_cell(cell_idx).cur_track = tracks{cell_idx};
        
        
    end
    
    
    
    return
    
    error('im not done yet.')
    
    % object_array is cell of linked objects
    % num_tracks is equal to number of actually linked objects
    
    linked_object_cell = cell(1,num_tracks);
    
    cur_linked_object_struct.track = [];
    cur_linked_object_struct.object_list = [];
    cur_linked_object_struct.frame_start = [];
    cur_linked_object_struct.frame_end = [];

    cur_linked_object_struct = repmat(cur_linked_object_struct,[1,num_tracks]);

    for frame_idx = 1:num_frames
        cur_objects = objects_cell{frame_idx};

        if isempty(cur_objects)
            continue
        end

        for object_idx = 1:numel(cur_objects)
            track_assignment_idx = find(object_ids(:,frame_idx) == object_idx);

            if numel(track_assignment_idx) > 1
                error('This object is assigned to more than one track')
            end

            if ~numel(track_assignment_idx)
                % this object is probably assigned to a track with less
                % than minimum_track_len frames

                break
        %                     error('No assignments to track for this object')
            end

            cur_linked_object_struct(linked_object_idx).object_list = ...
                [cur_linked_object_struct(linked_object_idx).object_list, cur_objects(object_idx)];
        end
    end

    cur_track = tracks{linked_object_idx};

    cur_linked_object_struct(linked_object_idx).track = cur_track;
    cur_linked_object_struct(linked_object_idx).frame_start = find(~isnan(sum(cur_track,1)),1,'first');
    cur_linked_object_struct(linked_object_idx).frame_end = find(~isnan(sum(cur_track,1)),1,'last');

    linked_object_cell = [linked_object_cell, {cur_linked_object_struct}];
    
%     for linked_object_idx = 1:num_tracks
%             
%         cur_linked_object_struct = [];
% 
%         cur_linked_object_struct.track = [];
%         cur_linked_object_struct.object_list = [];
%         cur_linked_object_struct.frame_start = [];
%         cur_linked_object_struct.frame_end = [];
% 
%         cur_linked_object_struct = repmat(cur_linked_object_struct,[1,num_tracks]);
% 
%         for frame_idx = 1:num_frames
%             cur_objects = objects_cell{frame_idx};
%             
%             if isempty(cur_objects)
%                 continue
%             end
% 
%             for object_idx = 1:numel(cur_objects)
%                 track_assignment_idx = find(object_ids(:,frame_idx) == object_idx);
%                 
%                 if numel(track_assignment_idx) > 1
%                     error('This object is assigned to more than one track')
%                 end
%                 
%                 if ~numel(track_assignment_idx)
%                     % this object is probably assigned to a track with less
%                     % than minimum_track_len frames
%                     
%                     break
% %                     error('No assignments to track for this object')
%                 end
%                 
%                 cur_linked_object_struct(linked_object_idx).object_list = ...
%                     [cur_linked_object_struct(linked_object_idx).object_list, cur_objects(object_idx)];
%             end
%         end
% 
%         cur_track = tracks{linked_object_idx};
% 
%         cur_linked_object_struct(linked_object_idx).track = cur_track;
%         cur_linked_object_struct(linked_object_idx).frame_start = find(~isnan(sum(cur_track,1)),1,'first');
%         cur_linked_object_struct(linked_object_idx).frame_end = find(~isnan(sum(cur_track,1)),1,'last');
% 
%         linked_object_cell = [linked_object_cell, {cur_linked_object_struct}];
%     end
end

