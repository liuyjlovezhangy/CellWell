function plotting(options)
    plopts = options.plot_options;
    
    im = importdata([options.filename '__analysis_results/input_movie.mat']);
    im = mat4D_to_gray(im);
    
    if plopts.well_segmentation_flag
        cw.plot.well_segmentation(well_segmentation_results_struct, im(:,:,:,options.bf_channel),...
            well_segmentation_results_struct.otsu_idcs, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
    end
    return
       
    % Plot well trajectories
    
    if plopts.well_tracks_flag
        movie_file = cw.plot.well_tracks(im, well_tracking_results_struct.im_shifted, well_segmentation_results_struct, well_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies'] );
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    % plot well fluorescence
    
    if plopts.well_by_well_flag
        movie_file = cw.plot.well_by_well(well_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);        
        
        if make_movies
            %copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell segmentation
        
    if plopts.cell_segmentation_flag
        movie_file = cw.plot.cell_segmentation_and_tracking(0,well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, [], signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell tracking
    
    if plopts.cell_tracking_flag
        movie_file = cw.plot.cell_segmentation_and_tracking('both',well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, cell_tracking_results_struct, signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end

    if plopts.cell_overlay_flag
        movie_file = cw.plot.cell_overlays(well_tracking_results_struct.wells, ...
            signal_detection_results_struct, cell_segmentation_results_struct,...
            cell_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
end




