function test_cell_tracking()

    %%% tracking parameters
    
    link_params.mode = 'conservative';
    link_params.searchrad = 10;
    link_params.gap_close = 0;
    
    link_params.min_track_len = 1;
    
    %%% load well information
    
    well_tracking_results_struct = importdata('code/test_cases/well_tracking_test_data.mat');
    signal_detection_results_struct = importdata('code/test_cases/noise_detection_test_data.mat');
    cell_segmentation_results_struct = importdata('code/test_cases/cell_segmentation_corrected_test_data.mat');
    
    %%% track cells
    
    num_frames = size(cell_segmentation_results_struct.detected_cell_props_water,2);
    
    link_params.total_num_frames = num_frames;
    
    cell_tracking_results_struct = track_cells( cell_segmentation_results_struct, link_params );
    
    save('code/test_cases/cell_tracking_test_data.mat','cell_tracking_results_struct');
    
    %%% draw results
    
    % draw linked object assignments
    
%     linked_objects_cells = cell_tracking_results_struct.linked_objects_cells;
%     
%     for obj_idx = 1:numel(linked_objects_cells)
%     
%         figure(13842)
%         clf
% 
%     end
    
    % plot everything else
    
    plot_cell_segmentation_and_tracking(1, well_tracking_results_struct.wells, cell_segmentation_results_struct, cell_tracking_results_struct, signal_detection_results_struct.is_noise_matrix, 1, '.');
%     plot_cell_overlays(well_tracking_results_struct.wells, signal_detection_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, 1, '.');
end

