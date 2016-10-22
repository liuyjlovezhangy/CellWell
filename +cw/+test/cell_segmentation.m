function test_cell_segmentation()
    % cell segmentation parameters

    cell_segmentation_opts.adaptive_thresh_scale = 1.5;
    cell_segmentation_opts.tracking_params.threshold_density = 0.1;
    cell_segmentation_opts.tracking_params.peak_stringency = 'low';
    cell_segmentation_opts.tracking_params.threshold_smoothing = 'off';
    
    cell_segmentation_opts.gauss_blur_on = 1;
    cell_segmentation_opts.gauss_blur_sigma = 3;
    cell_segmentation_opts.gauss_blur_hsize = [5 5];

    cell_segmentation_opts.min_cell_area = 15;
    
    cell_segmentation_opts.watershedding = [1 1 1];
    cell_segmentation_opts.close_radius = [1 1 4];
    
    % load well information
    
    well_tracking_results_struct = importdata('code/test_cases/well_tracking_test_data.mat');
    signal_detection_results_struct = importdata('code/test_cases/noise_detection_test_data.mat');
    
    % segment cells
    
    cell_segmentation_results_struct = segment_cells( well_tracking_results_struct, signal_detection_results_struct, cell_segmentation_opts );
    
    save('code/test_cases/cell_segmentation_test_data.mat','cell_segmentation_results_struct');
    
    % draw results
    
%     plot_cell_segmentation_and_tracking(0, well_tracking_results_struct.wells, cell_segmentation_results_struct, [], signal_detection_results_struct.is_noise_matrix, 1, '.');
end

