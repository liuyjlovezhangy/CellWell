function nuclei_cfg

    options.filename = 'movies/Farid_test-02(1)_NEW.tif';
    
    movie_name = strrep(options.filename,'/','__');
    
    %%% RUNNING MODE
    
    options.DO_PROCESSING = 0;
    options.DO_ANALYSIS = 1;
    options.DO_PLOTTING = 0;
    
    options.NUKE_IT = 0; %%% DESTROY ALL OLD OUTPUT DATA/MOVIES INCL OLD RUNS
    options.NUKE_WARN = 0; %%% Pops up a box every time when run with NUKE_IT=1 to ask
    
    options.ask_me = 1; %%% Ask user if the results look good after registration, segmentation, etc.
    
    options.first_frame = 1;
    options.last_frame = 40; %%% Only process through <n> frames
    
    %%% Info related to movie
    
    options.pixel_size = 1/0.75; % microns per pixel
    options.time_step = 3*60; % Number of seconds per frame
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CHANGE THESE PER MOVIE SETUP
    
    %%% ONE SIGNAL CHANNEL, TWO CELL CHANNELS, 75 um WELLS
    
    % Movie channel setup
    
    options.signal_channels = [];
    options.cell_channels = [1 2];
    options.nuclei_channel = [];
    options.tracking_channels = [1 2];
    options.bf_channel = 4;
    
    options.channel_labels = {'T cell','Target','Nuclei','Brightfield'};
    
    % Plate / well setup
    
    options.well_counts = [6,6]; % #rows x #cols
    options.well_width = 69;
    options.well_spacing_width = 4;
    options.well_height = 72;
    options.well_spacing_height = 2;
    
    % cell information
    
    options.interaction_tolerance = 6; % how many pixels apart cell edges can be to consider them interacting
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% BASIC processing options for segmenting wells, tracking cells, etc.
    %%%
    
    %%% What part of the processing would you like to start from?
    % Uncomment the line you are interested in (and comment out others)
    
%      processing_options.start_at = 0;
%     processing_options.start_at = 'well segmentation';
%     processing_options.start_at = 'well tracking';
%     processing_options.start_at = 'noise detection';
    processing_options.start_at = 'cell segmentation';
%     processing_options.start_at = 'cell tracking';
%     processing_options.start_at = 'next';
    
    %%% Movie modifications
    
    processing_options.register = 0;
    processing_options.rotate = 1;

    %%% Well segmentation mode
    
    processing_options.wseg_debug = 0;
    processing_options.wseg_mode = 'template';
    processing_options.wseg_template = 'well_templates/wt_default.tif';
    processing_options.wseg_left_offset = 5;
    processing_options.wseg_top_offset = 5;
    
    %%% Cell segmentation mode
    
    processing_options.cseg_debug = 0;
    processing_options.cseg_mode = 'radial';
    
    %%% Cell object cleanup
    %%% What is the minimum area of a cell object, etc.
    
    processing_options.cseg_min_area = 0;

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
    
    plot_options.well_segmentation = 0;
    plot_options.well_tracks = 0;
    plot_options.well_by_well = 0;
    plot_options.noise_detection = 0;
    plot_options.cell_segmentation = 0; % BROKEN FOR RADIAL
    plot_options.cell_tracking = 0; % BROKEN FOR RADIAL cell tracking also plots cell segmentation.
    plot_options.cell_overlay = 1;
    
    plot_options.make_movies = 1;
    
    %%% Define what wells are "interesting" for plotting
    
    plot_options.min_cells = [0 0 0]; % how many cells necessary in channel 1, 2, ...
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% ADVANCED Options for segmenting wells, tracking cells, etc.
    %%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Well segmentation parameters
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Noise detection parameters
    
    processing_options.ndec_thresh_mean = 2;
    processing_options.ndec_thresh_stdev = 0.5;
    processing_options.ndec_do_wiener = 1;
    processing_options.ndec_do_blur = 1;
    processing_options.ndec_blur_sigma = 10;
    processing_options.ndec_blur_hsize = [10 10];
    processing_options.ndec_entropy_nsize = ones(9);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell segmentation parameters
    
    %%% GENERAL THRESHOLDING OPTIONS
    
    processing_options.cseg_adaptive_thresh_scale = 0.6;
    processing_options.cseg_tracking_params.threshold_density = 0.1;
    processing_options.cseg_tracking_params.peak_stringency = 'low';
    processing_options.cseg_tracking_params.threshold_smoothing = 'off';
    
    processing_options.cseg_min_cell_area = 20;   
    
    %%% Cell tracking parameters
    
    processing_options.ctrk_mode = 'conservative';
    processing_options.ctrk_searchrad = 10;
    processing_options.ctrk_gap_close = 0;
    
    processing_options.ctrk_min_track_len = 3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Compile options and run
    %%%
    
    options.processing_options = processing_options;
    options.analysis_options = analysis_options;
    options.plot_options = plot_options;
    
    save([movie_name '_cfg.mat'],'options')

end