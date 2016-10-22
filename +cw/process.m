function process(options)

    propts = options.processing_options;

    disp('PROCESSING MOVIE...')

    multiWaitbar('CloseAll');
    
    full_filename = options.filename;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Load movie 

    if ~propts.start_at || ~exist([full_filename '__analysis_results/input_movie.mat'],'file')
        disp('Loading movie from image file...')
        
        im = zloadim(full_filename,1);
        
        save([full_filename '__analysis_results/input_movie.mat'],'im')
        
    else
        disp('Loading movie from mat...')
        
        im = importdata([full_filename '__analysis_results/input_movie.mat']);
    end
    
    im = mat4D_to_gray(im);
    
%     link_params.total_num_frames = size(im,3);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set up data structures
    
    well_segmentation_results_struct = [];
    frame = [];
    
    frame.mask = [];
    frame = repmat(frame,[1,size(im,3)]);
    
    well_segmentation_results_struct.frame = frame;
    well_segmentation_results_struct.well_tracks = [];
    well_segmentation_results_struct.well_ids = [];
    well_segmentation_results_struct.extra_border_x = propts.wseg_extra_border_x;
    well_segmentation_results_struct.extra_border_y = propts.wseg_extra_border_y;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform segmentation on each BF frame
    
    if propts.start_at <= 1 || ~exist([full_filename '__analysis_results/well_segmentation.mat'],'file')
    
        well_segmentation_results_struct = cw.process.segment_wells( im(:,:,:,propts.bf_channel), propts );
        
        disp('Saving well segmentation results...')
        
        save([full_filename '__analysis_results/well_segmentation.mat'],'well_segmentation_results_struct')
        %copyfile([full_filename '__analysis_results/well_segmentation.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/well_segmentation.mat']);
    else
        disp('Loading previous well segmentation...')
        
        well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    end
    
    return
    
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
end




