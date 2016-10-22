function test_cfg
    
    options.filename = 'movies/Timelapse01.ome.tiff';

    %%% What part of the analysis would you like to start from?
    % Uncomment the line you are interested in (and comment out others)
    
    options.start_at = 0;
%     options.start_at = 'well segmentation';
%     options.start_at = 'well tracking';
%     options.start_at = 'noise detection';
%     options.start_at = 'cell segmentation';
%     options.start_at = 'cell tracking';
%     options.start_at = 'analysis';
    
    %%% Info related to movie
    
    options.pixel_size = 1/0.75; % microns per pixel
    options.time_step = 5*60; % Number of seconds per frame
    
    %%% Movie channel setup
    
    options.signal_channel = 1;
    options.cell_channels = [2 3];
    options.tracking_channels = [2 3];
    options.bf_channel = 4;

    %%% PLOTTING
    
    options.plot_well_segmentation = 0;
    options.plot_well_tracks = 0;
    options.plot_well_by_well = 0;
    options.plot_noise_detection = 0;
    options.plot_cell_segmentation = 0;
    options.plot_cell_tracking = 0;
    options.plot_cell_overlay = 0;
    
    options.make_movies = 1;
    
    %%% MISC
    
    options.correct_rotation = 1; %%% WARNING: NOT RIGOROUSLY TESTED
    
    options
    
%     cw.process(options);
end