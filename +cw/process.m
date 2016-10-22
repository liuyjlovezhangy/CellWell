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
    %%% Perform segmentation on each BF frame
    
    if propts.start_at <= 1 || ~exist([full_filename '__analysis_results/well_segmentation.mat'],'file')
    
        disp('Performing well segmentation...')
        
        well_segmentation_results_struct = cw.process.segment_wells( im(:,:,:,propts.bf_channel), propts );
        
        disp('Saving well segmentation results...')
        
        save([full_filename '__analysis_results/well_segmentation.mat'],'well_segmentation_results_struct')
        %copyfile([full_filename '__analysis_results/well_segmentation.mat'],[full_filename '__analysis_results/runs/' analysis_timestamp,'/well_segmentation.mat']);
    else
        disp('Loading previous well segmentation...')
        
        well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform well tracking to link wells from frame to frame
    
    if propts.start_at <= 2 || ~exist([full_filename '__analysis_results/well_tracking.mat'],'file')
    
        disp('Performing well tracking...')
        
        well_tracking_results_struct = cw.process.track_wells( im, well_segmentation_results_struct, propts );
    
        disp('Saving well tracking results...')
        
        save([full_filename '__analysis_results/well_tracking.mat'],'well_tracking_results_struct','-v7.3')

    else
        disp('Loading previous well tracking...')
        
        well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
    end
    
    return
    
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




