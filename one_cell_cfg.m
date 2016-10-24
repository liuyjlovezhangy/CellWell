function one_cell_cfg

    options.filename = 'movies/CARTonly_75um_before_leak.tif';
    
    %%% RUNNING MODE
    
    options.DO_PROCESSING = 0;
    options.DO_ANALYSIS = 0;
    options.DO_PLOTTING = 1;
    
    options.NUKE_IT = 0; %%% DESTROY ALL OLD OUTPUT DATA/MOVIES INCL OLD RUNS
    options.NUKE_WARN = 0; %%% Pops up a box every time when run with NUKE_IT=1 to ask
    
    %%% Info related to movie
    
    options.pixel_size = 1/0.75; % microns per pixel
    options.time_step = 5*60; % Number of seconds per frame
    
    % Movie channel setup
    
    options.signal_channel = [];
    options.cell_channels = 1;
    options.tracking_channels = 1;
    options.bf_channel = 2;
    
    options.channel_names = {'(CAR)T cell','Brightfield'};
    
    % Plate / well setup
    
    % 
    
    options.well_counts = [4,4];
    options.well_width = 97;
    options.well_spacing_width = 14;
    options.well_height = 98;
    options.well_spacing_height = 13;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Options for segmenting wells, tracking cells, etc.
    %%%

    %%% What part of the processing would you like to start from?
    % Uncomment the line you are interested in (and comment out others)
    
    processing_options.start_at = 0;
%     processing_options.start_at = 'well segmentation';
%     processing_options.start_at = 'well tracking';
%     processing_options.start_at = 'noise detection';
%     processing_options.start_at = 'cell segmentation';
%     processing_options.start_at = 'cell tracking';
%     processing_options.start_at = 'next';
    
    %%% Well segmentation parameters
    
    % THESE OPTIONS ARE FOR THE ORIGINAL (FAST BUT NOT ROBUST) SEGMENTATION
    
    processing_options.wseg_mode = 'otsu';
    
    processing_options.wseg_otsu_close_size = 10;
    processing_options.wseg_otsu_min_size = 2000;
    processing_options.wseg_otsu_num_otsu_levels = 2;
    processing_options.wseg_otsu_extra_border_x = 3;
    processing_options.wseg_otsu_extra_border_y = 5;
    
    % THE OPTIONS ARE FOR NEW (SLOW BUT ROBUST) SEGMENTATION
    
    processing_options.wseg_mode = 'edge';
    
    processing_options.wseg_edge_
    
    %%% Noise detection parameters
    
    processing_options.ndec_thresh_mean = 1;
    processing_options.ndec_thresh_stdev = 0.5;
    processing_options.ndec_do_wiener = 1;
    processing_options.ndec_do_blur = 1;
    processing_options.ndec_blur_sigma = 10;
    processing_options.ndec_blur_hsize = [10 10];
    processing_options.ndec_entropy_nsize = ones(9);
    
    %%% Cell segmentation parameters
    
    processing_options.cseg_adaptive_thresh_scale = 1.5;
    processing_options.cseg_tracking_params.threshold_density = 0.1;
    processing_options.cseg_tracking_params.peak_stringency = 'low';
    processing_options.cseg_tracking_params.threshold_smoothing = 'off';

    processing_options.cseg_min_cell_area = 10;
    
    processing_options.cseg_watershedding = [1 1 1];
    processing_options.cseg_close_radius = [1 1 4];
    
    %%% Cell tracking parameters
    
    processing_options.ctrk_mode = 'conservative';
    processing_options.ctrk_searchrad = 10;
    processing_options.ctrk_gap_close = 0;
    
    processing_options.ctrk_min_track_len = 3;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Options for data analysis post-processing
    %%%
    
    analysis_options = [];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Options for plotting
    %%%
    
    plot_options.well_segmentation = 1;
    plot_options.well_tracks = 0;
    plot_options.well_by_well = 0;
    plot_options.noise_detection = 0;
    plot_options.cell_segmentation = 0;
    plot_options.cell_tracking = 0; % cell tracking also plots cell segmentation.
    plot_options.cell_overlay = 0;
    
    plot_options.make_movies = 1;

    %%% Compile options and run
    
    options.processing_options = processing_options;
    options.analysis_options = analysis_options;
    options.plot_options = plot_options;
    
    cw.main(options);
end