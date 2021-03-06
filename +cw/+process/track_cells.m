function cell_tracking_results_struct = track_cells( cell_segmentation_results_struct, options )

    propts = options.processing_options;

    num_wells = numel(cell_segmentation_results_struct.cell_masks);
    num_channels = size(cell_segmentation_results_struct.detected_cell_props,3);
    num_frames = size(cell_segmentation_results_struct.detected_cell_props,2);
    
    cell_tracking_results_struct = [];
    cell_tracking_results_struct.cell_tracks = cell(num_wells,num_channels);
    cell_tracking_results_struct.cell_ids = cell(num_wells,num_channels);
    cell_tracking_results_struct.linked_object_cells = cell(num_wells,num_channels);

    link_params.mode = propts.ctrk_mode;
    link_params.searchrad = propts.ctrk_searchrad;
    link_params.gap_close = propts.ctrk_gap_close;
    link_params.total_num_frames = num_frames;
    
    % build a list of all cell locations throughout time

    for well_idx = 1:num_wells
        for channel_idx = options.tracking_channels
            objects_cell = cell(1,num_frames);

            localization_array = [];

            for frame_idx = 1:num_frames

                cur_objects = cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,channel_idx};

                objects_cell{frame_idx} = cur_objects;

                for obj_idx = 1:numel(cur_objects)
                    localization_array = [localization_array; [cur_objects(obj_idx).Centroid, frame_idx, obj_idx]];
                end
            end

            if ~isempty(localization_array)
                %%% link the cell positions frame-to-frame

                [cell_tracks,cell_ids] = cw.process.simple_tracking(localization_array, link_params );

                %%% match the regionprops objects with their tracks
                
                linked_object_cell = cw.process.assign_tracks_to_objects( cell_tracks, cell_ids, objects_cell );
                
                cell_tracking_results_struct.cell_tracks{well_idx,channel_idx} = cell_tracks;
                cell_tracking_results_struct.linked_object_cells{well_idx,channel_idx} = linked_object_cell;
                cell_tracking_results_struct.cell_ids{well_idx,channel_idx} = cell_ids;
            end
        end
        
    end
    
    cell_tracking_results_struct.link_params = link_params;
end

