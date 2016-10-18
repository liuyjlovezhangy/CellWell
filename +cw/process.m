function process(full_filename, start_fresh_flag)

    multiWaitbar('CloseAll');

    if ~exist(full_filename,'file')
        msgbox('ERROR: Could not find supplied movie file')
    end

    if ~exist('start_fresh_flag','var')
        start_fresh_flag = 0;
    end

    % Plotting options
    
    plot_well_segmentation_flag = 1;
    plot_well_tracks_flag = 1;
    plot_well_by_well_flag = 1;
    plot_noise_detection_flag = 1;
    plot_cell_segmentation_flag = 1;
    plot_cell_tracking_flag = 1;
    plot_cell_overlay_flag = 1;
    plot_statistics_flag = 1;
    
    make_movies = 1;
    
    % Imaging options
    %%%%%%%%%%%%%%% THIS SHOULD BE USER DEFINED PER MOVIE
    signal_channel = 1;
    bf_channel = 4; 
    
    % Well segmentation options
    
    num_otsu_levels = 2;

    % How much bigger are the wells than the segmentation mask finds?
    
    extra_border_x = 3;
    extra_border_y = 5;

    % noise detection parameters
    
    detection_opts.thresh_mean = 1;
    detection_opts.thresh_stdev = 0.5;
    detection_opts.do_wiener = 1;
    detection_opts.do_blur = 1;
    detection_opts.blur_sigma = 10;
    detection_opts.blur_hsize = [10 10];
    detection_opts.entropy_nsize = ones(9);
    
    % options to decide if a well was segmented correctly
    
    [~,o] = cw.analyze.well_criterion();
    
    % cell segmentation
    
    cell_segmentation_opts.adaptive_thresh_scale = 1.5;
    cell_segmentation_opts.tracking_params.threshold_density = 0.1;
    cell_segmentation_opts.tracking_params.peak_stringency = 'low';
    cell_segmentation_opts.tracking_params.threshold_smoothing = 'off';

    cell_segmentation_opts.min_cell_area = 10;
    
    cell_segmentation_opts.watershedding = [1 1 1];
    cell_segmentation_opts.close_radius = [1 1 4];
    
    % cell tracking
    
    link_params.mode = 'conservative';
    link_params.searchrad = 10;
    link_params.gap_close = 0;
    
    link_params.min_track_len = 3;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Directory setup and detection of progress
    
    [path_name,filename,ext] = fileparts(full_filename);
    
    mkdir([full_filename '__analysis_results'])
    mkdir([full_filename '__analysis_results/movies'])
    mkdir([full_filename '__analysis_results/runs'])
    
%     analysis_timestamp =
%     strrep(datestr(datetime('now','TimeZone','local','Format','d-MMM-y HH:mm:ss')),':','-'); % only supported by r2014b and above
    analysis_timestamp = strrep(datestr(now),':','-');    
    mkdir([full_filename '__analysis_results/runs/' analysis_timestamp])
    mkdir([full_filename '__analysis_results/runs/' analysis_timestamp '/movies'])
    
    if start_fresh_flag
        delete([full_filename '__analysis_results/*.avi'])
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load movie 

    if ~exist([full_filename '__analysis_results/movie.mat'],'file') || start_fresh_flag
        disp('Loading movie from image file...')
        
        im = zloadim(full_filename,1);
        
        % Make it so the bf channel is the last one
        % Signal channel should be blue

        channel_order = 1:size(im,4);
        channel_order((channel_order == bf_channel) | (channel_order == signal_channel)) = [];
        channel_order = [channel_order, signal_channel, bf_channel];

        channel_order(channel_order == bf_channel) = [];
        channel_order = [channel_order, bf_channel];
        im = im(:,:,:,channel_order);
        
        save([full_filename '__analysis_results/movie.mat'],'im')
        
    else
        disp('Loading movie from mat...')
        
        im = importdata([full_filename '__analysis_results/movie.mat']);
    end
    
    im = mat4D_to_gray(im);
    
    link_params.total_num_frames = size(im,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set up data structures
    
    well_segmentation_results_struct = [];
    frame = [];
    
    frame.mask = [];
    frame = repmat(frame,[1,size(im,3)]);
    
    well_segmentation_results_struct.frame = frame;
    well_segmentation_results_struct.well_tracks = [];
    well_segmentation_results_struct.well_ids = [];
    well_segmentation_results_struct.extra_border_x = extra_border_x;
    well_segmentation_results_struct.extra_border_y = extra_border_y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform segmentation on each BF frame
    
    if ~exist([full_filename '__analysis_results/well_segmentation.mat'],'file') || start_fresh_flag
    
        disp('Performing well segmentation...')
        
        h=waitbar(0,'Segmenting movie into wells...') ;
    
        for frame_idx = 1:size(im,3)

            im_bf_slice = im(:,:,frame_idx,end);

            %%% Find otsu threshold levels and binarize

            [IDX,sep,well_segmentation_results_struct.otsu_idcs] = otsu(im_bf_slice,num_otsu_levels);

            im_mask_thresh = IDX == num_otsu_levels;

            %%% Image morphological operations to clean up

            % remove objects on border

            im_mask_noborder = imclearborder(im_mask_thresh);

            % link together objects

            se = strel('disk',10);
            im_mask_connect = imclose(im_mask_noborder,se);

            % delete border objects again

            im_mask_noborder2 = imclearborder(im_mask_connect);

            % Clear very small objects
            im_mask_cleared = bwareaopen(im_mask_noborder2, 0.9*o.min_well_area);

            im_mask_final = im_mask_cleared;
            im_mask_final_corrected = im_mask_final;

            im_seg = im_bf_slice .* im_mask_final;

            %%% Object detection

            L = bwlabeln(im_seg);
            objects = regionprops(L,im_seg,'all');

            good_objects = [];

            %%% Cull bad objects

            for obj_idx = 1:numel(objects)
                if cw.analyze.well_criterion(objects(obj_idx))
                    good_objects = [good_objects, objects(obj_idx)];
                    
                    extrema = objects(obj_idx).Extrema;
                
                    xmin = floor(min(extrema(:,1)) - extra_border_x);
                    xmax = ceil(max(extrema(:,1)) + extra_border_x);
                    ymin = floor(min(extrema(:,2)) - extra_border_y);
                    ymax = ceil(max(extrema(:,2)) + extra_border_y);
                    
                    im_mask_final_corrected(ymin:ymax,xmin:xmax) = 1;
                end
            end

            well_segmentation_results_struct.frame(frame_idx).all_objects = objects;
            well_segmentation_results_struct.frame(frame_idx).good_objects = good_objects;
            well_segmentation_results_struct.frame(frame_idx).mask = im_mask_final;
            well_segmentation_results_struct.frame(frame_idx).im_mask_thresh = im_mask_thresh;
            well_segmentation_results_struct.frame(frame_idx).im_mask_final_corrected = im_mask_final_corrected;

            % create new mask and segmented image based on programmatic
            % well borders
            
            waitbar(frame_idx / size(im,3))
            drawnow
        end

        close(h);
        drawnow
        
        disp('Saving well segmentation results...')
        
        save([full_filename '__analysis_results/well_segmentation.mat'],'well_segmentation_results_struct')
        %copyfile([full_filename '__analysis_results/well_segmentation.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/well_segmentation.mat']);
    else
        disp('Loading previous well segmentation...')
        
        well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    end
    
    if plot_well_segmentation_flag
        movie_file = cw.plot.well_segmentation(well_segmentation_results_struct, im(:,:,:,end), well_segmentation_results_struct.otsu_idcs, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform well tracking to link wells from frame to frame
    
    if ~exist([full_filename '__analysis_results/well_tracking.mat'],'file') || start_fresh_flag
    
        disp('Performing well tracking...')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Initial tracking

        % build a list of all well locations throughout time

        localization_array = [];

        for frame_idx = 1:size(im,3)
            for well_idx = 1:numel(well_segmentation_results_struct.frame(frame_idx).good_objects)
                objects = well_segmentation_results_struct.frame(frame_idx).good_objects(well_idx);

                localization_array = [localization_array; [objects.Centroid, frame_idx, well_idx]];

            end
        end

        % link the well positions frame-to-frame

        linkparams.mode = 'conservative';
        linkparams.searchrad = 20;
        linkparams.gap_close = 1;

        [well_tracks,well_ids] = cw.analyze.simple_tracking( localization_array, linkparams );

        well_tracking_results_struct.well_tracks = well_tracks;
        well_tracking_results_struct.well_ids = well_ids;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Correct for stage drift

        % Get average position trajectory

        well_tracks_mat = tracks_to_matrix(well_tracks);
        well_tracks_mat_shifted = zeros(size(well_tracks_mat));

        starting_average_position = squeeze(mean(well_tracks_mat,1));

        im_shifted = zeros(size(im));

        for frame_idx = 1:size(im,3)
            
            im_slice_cur = squeeze(im(:,:,frame_idx,:));

            % calculate shift

            shift = starting_average_position(frame_idx,:) - starting_average_position(1,:);
            im_shift = fliplr(shift);

            % shift image
            im_shifted(:,:,frame_idx,:) = shift_image(im_slice_cur,-im_shift(1),-im_shift(2));

            % correct tracks

            well_tracks_mat_shifted(:,frame_idx,:) = squeeze(well_tracks_mat(:,frame_idx,:)) - repmat(shift,[size(well_tracks_mat,1), 1]);
        end

        well_tracks_shifted = {};

        for track_idx = 1:size(well_tracks_mat_shifted,1)
            well_tracks_shifted = [well_tracks_shifted, {squeeze(well_tracks_mat_shifted(track_idx,:,:))'}];
        end

        well_tracking_results_struct.well_tracks_shifted = well_tracks_shifted;

        well_tracking_results_struct.well_tracks_mat = well_tracks_mat;
        well_tracking_results_struct.well_tracks_mat_shifted = well_tracks_mat_shifted;

        well_tracking_results_struct.im_shifted = im_shifted;

        % Add a wells entry to the structure where r_s.wells(well_idx) is
        % itself a structure containing the actual well's trajectory and
        % others

        total_wells = max(well_ids(:));

        wells.track = [];

        % measure the average x,y length of the wells in the first frame

        well_half_widths = [];
        well_half_heights = [];

        for obj_idx = 1:numel(well_segmentation_results_struct.frame(1).good_objects)
            props = well_segmentation_results_struct.frame(1).good_objects(obj_idx);

            well_half_widths = [well_half_widths, ceil((max(props.Extrema(:,1)) - min(props.Extrema(:,1))) ./ 2)];
            well_half_heights = [well_half_heights, ceil((max(props.Extrema(:,2)) - min(props.Extrema(:,2))) ./ 2)];
        end

        well_half_widths = mean(well_half_widths);
        well_half_heights = mean(well_half_heights);
        
        new_mask = zeros(size(im,1),size(im,2));

        wells = repmat(wells,[1,total_wells]);

        for well_idx = 1:total_wells
            wells(well_idx).track = well_tracks{well_idx};
            wells(well_idx).track_shifted = well_tracks_shifted{well_idx};

            wells(well_idx).left_boundary = mean(well_tracks_shifted{well_idx}(1,:)) - well_half_widths - extra_border_x;
            wells(well_idx).right_boundary = mean(well_tracks_shifted{well_idx}(1,:)) + well_half_widths + extra_border_x;
            wells(well_idx).bottom_boundary = mean(well_tracks_shifted{well_idx}(2,:)) - well_half_heights - extra_border_y;
            wells(well_idx).top_boundary = mean(well_tracks_shifted{well_idx}(2,:)) + well_half_heights + extra_border_y;

            % make new mask now based on well tracking
            
            i_idcs = floor(wells(well_idx).bottom_boundary:ceil(wells(well_idx).top_boundary));
            j_idcs = floor(wells(well_idx).left_boundary):ceil(wells(well_idx).right_boundary);

            new_mask(i_idcs,j_idcs) = 1;
            
%             wells(well_idx).im_well = gray_to_uint16(im_shifted(i_idcs,j_idcs,:,:));
            wells(well_idx).im_well = im_shifted(i_idcs,j_idcs,:,:);
            
        end

        well_tracking_results_struct.wells = wells;
        well_tracking_results_struct.mask = new_mask;
    
        disp('Saving well tracking results...')
        
%         well_tracking_results_struct.im_shifted = gray_to_uint16(well_tracking_results_struct.im_shifted);
        
        save([full_filename '__analysis_results/well_tracking.mat'],'well_tracking_results_struct','-v7.3')
        %copyfile([full_filename '__analysis_results/well_tracking.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/well_tracking.mat']);

    else
        disp('Loading previous well tracking...')
        
        well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
    end
    
%     for well_idx = 1:numel(well_tracking_results_struct.wells)
%         well_tracking_results_struct.wells(well_idx).im_well = mat4D_to_gray(well_tracking_results_struct.wells(well_idx).im_well);
%     end
    
%     well_tracking_results_struct.im_shifted = mat4D_to_gray(well_tracking_results_struct.im_shifted);
    
    % Plot well trajectories
    
    if plot_well_tracks_flag
        movie_file = cw.plot.well_tracks(im, well_tracking_results_struct.im_shifted, well_segmentation_results_struct, well_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies'] );
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    % plot well fluorescence
    
    if plot_well_by_well_flag
        movie_file = cw.plot.well_by_well(well_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);        
        
        if make_movies
            %copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Well-by-well detection of signal presence
    
    if ~exist([full_filename '__analysis_results/noise_detection.mat'],'file') || start_fresh_flag
        
        disp('Detecting wells and frames with signal in each channel...')
        
        h = waitbar(0,'Detecting wells and frames that contain signal [cells]...');

        num_wells = numel(well_tracking_results_struct.wells);

        is_noise_matrix = zeros(num_wells,size(im,4)-1,size(im,3)); % well, channel, frame
        detection_images = cell(1,num_wells);

        for well_idx = 1:numel(well_tracking_results_struct.wells)
            im_well = well_tracking_results_struct.wells(well_idx).im_well;

            % entropy calc chops off pixels on img border due to artifacting
            [~,~,sizing_im] = cw.analyze.detect_noise( im_well(:,:,1,1),detection_opts );
            cur_detection_image = zeros(size(sizing_im,1), size(sizing_im,2), size(im_well,3), size(im,4)-1); 

            for channel_idx = 1:size(im,4)-1 

                [temp_noise_matrix,~,channel_im_filtered] = cw.analyze.detect_noise( im_well(:,:,:,channel_idx), detection_opts );

                is_noise_matrix(well_idx,channel_idx,:) = temp_noise_matrix;
                cur_detection_image(:,:,:,channel_idx) = channel_im_filtered;

            end

            detection_images{well_idx} = cur_detection_image;

            waitbar(well_idx / numel(well_tracking_results_struct.wells));
        end

        signal_detection_results_struct.detection_images = detection_images;
        signal_detection_results_struct.is_noise_matrix = is_noise_matrix;
        signal_detection_results_struct.detection_opts = detection_opts;

        close(h);
        drawnow
        
        disp('Saving signal detection results....')
        
        save([full_filename '__analysis_results/noise_detection.mat'],'signal_detection_results_struct')
        %copyfile([full_filename '__analysis_results/noise_detection.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/noise_detection.mat']);
        
    else
        disp('Loading previous signal detection data...')
        
        signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
    end
    
    if plot_noise_detection_flag
        movie_file = cw.plot.signal_detection( well_tracking_results_struct.wells,signal_detection_results_struct.detection_images,...
            signal_detection_results_struct.is_noise_matrix,signal_detection_results_struct.detection_opts.thresh_mean,signal_detection_results_struct.detection_opts.thresh_stdev,make_movies,[full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell segmentation
        
    if ~exist([full_filename '__analysis_results/cell_segmentation.mat'],'file') || start_fresh_flag
        
        disp('Segmenting cells in each well...')

        cell_segmentation_results_struct = cw.analyze.segment_cells( well_tracking_results_struct, signal_detection_results_struct, cell_segmentation_opts );
        
        disp('Cleaning up segmentation....')
        cell_segmentation_results_struct = cw.analyze.clean_up_cell_segmentation( cell_segmentation_results_struct, [] );
       
        disp('Saving cell segmentation results....')
        
        save([full_filename '__analysis_results/cell_segmentation.mat'],'cell_segmentation_results_struct')
        %copyfile([full_filename '__analysis_results/cell_segmentation.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/cell_segmentation.mat'])
    else
        disp('Loading previous cell segmentation data...')
        
        cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        
    end
    
    if plot_cell_segmentation_flag
        movie_file = cw.plot.cell_segmentation_and_tracking(0,well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, [], signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell tracking
    
    if ~exist([full_filename '__analysis_results/cell_tracking.mat'],'file') || start_fresh_flag
        
        disp('Tracking cells in each well...')
        
        cell_tracking_results_struct = cw.analyze.track_cells( cell_segmentation_results_struct, link_params );
        
        % save
        
        disp('Saving cell tracking results....')
        
        save([full_filename '__analysis_results/cell_tracking.mat'],'cell_tracking_results_struct')
        %copyfile([full_filename '__analysis_results/cell_tracking.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/cell_tracking.mat'])
    else
        disp('Loading previous cell tracking...')
        
        cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    end
    
    if plot_cell_tracking_flag
        movie_file = cw.plot.cell_segmentation_and_tracking('both',well_tracking_results_struct.wells, ...
            cell_segmentation_results_struct, cell_tracking_results_struct, signal_detection_results_struct.is_noise_matrix, make_movies,...
            [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end

    if plot_cell_overlay_flag
        movie_file = cw.plot.cell_overlays(well_tracking_results_struct.wells, ...
            signal_detection_results_struct, cell_segmentation_results_struct,...
            cell_tracking_results_struct, make_movies, [full_filename '__analysis_results/runs/' analysis_timestamp '/movies']);
        
        if make_movies
            copyfile(movie_file,[full_filename '__analysis_results/movies/']);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Data analysis
    
    cw.analyze.final_analysis(well_tracking_results_struct,signal_detection_results_struct,cell_segmentation_results_struct,cell_tracking_results_struct);
end




