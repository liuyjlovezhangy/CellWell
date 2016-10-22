function plotting(options)
    plopts = options.plot_options;
    make_movies = plopts.make_movies;
    
    full_filename = options.filename;
    
    %%% Load movie
    
    im = importdata([full_filename '__analysis_results/input_movie.mat']);
    im = mat4D_to_gray(im);
    
    well_segmentation_results_struct = [];
    well_tracking_results_struct = [];
    signal_detection_results_struct = [];
    cell_segmentation_results_struct = [];
    cell_tracking_results_struct = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Well segmentation-related movies
    
    if plopts.well_segmentation
        if isempty(well_segmentation_results_struct)
            well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
        end
        
        cw.plot.well_segmentation(well_segmentation_results_struct, im(:,:,:,options.bf_channel),...
            make_movies, [full_filename '__analysis_results/movies']);
    end
    
    % Plot well trajectories
    
    if plopts.well_tracks
        if isempty(well_segmentation_results_struct)
            well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
        end
        
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        cw.plot.well_tracks(im, well_tracking_results_struct.im_shifted, well_segmentation_results_struct, ...
            well_tracking_results_struct, make_movies, [full_filename '__analysis_results/movies'] );

    end
    
    % plot well fluorescence
    
    if plopts.well_by_well
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        cw.plot.well_by_well(well_tracking_results_struct, make_movies,...
            [full_filename '__analysis_results/movies']);        
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Single well cell movies
        
    % signal
    
    if plopts.noise_detection
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        if isempty(signal_detection_results_struct)
           signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
        end
        
        cw.plot.signal_detection( well_tracking_results_struct.wells,signal_detection_results_struct.detection_images,...
                signal_detection_results_struct.is_noise_matrix,signal_detection_results_struct.detection_opts.ndec_thresh_mean,...
                signal_detection_results_struct.detection_opts.ndec_thresh_stdev,make_movies,[full_filename '__analysis_results/movies']);
    end
    
    % segmentation
    
    if plopts.cell_segmentation
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        if isempty(cell_segmentation_results_struct)
           cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        end
        
        if isempty(signal_detection_results_struct)
           signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
        end
        
        cw.plot.cell_segmentation_and_tracking(0,well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, [], signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/movies']);
    end

    % tracking
    
    if plopts.cell_tracking
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        if isempty(cell_segmentation_results_struct)
           cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        end
        
        if isempty(cell_tracking_results_struct)
           cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
        end
        
        if isempty(signal_detection_results_struct)
           signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
        end
        
        cw.plot.cell_segmentation_and_tracking('both',well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, cell_tracking_results_struct, signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/movies']);
    end

    % overlay
    
    if plopts.cell_overlay
        if isempty(well_tracking_results_struct)
           well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
        end
        
        if isempty(cell_segmentation_results_struct)
           cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        end
        
        if isempty(cell_tracking_results_struct)
           cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
        end
        
        if isempty(signal_detection_results_struct)
           signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
        end
        
        cw.plot.cell_overlays(well_tracking_results_struct.wells, ...
            signal_detection_results_struct, cell_segmentation_results_struct,...
            cell_tracking_results_struct, make_movies, [full_filename '__analysis_results/movies']);
    end
end




