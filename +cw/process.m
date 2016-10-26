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
        
        % Perform registration and rotation if desired
        
        if options.processing_options.register
            im_reg = cw.process.register(im,options.bf_channel);
            
            if options.ask_me

                im_combination = [im(:,:,:,options.bf_channel), im_reg(:,:,:,options.bf_channel)];

                answer = cw.plot.confirm_results(im_combination(:,:,:,options.bf_channel),'Registration results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end
            
            im = im_reg;
        end
        
        if options.processing_options.rotate
            im_rotate = cw.process.rotate(im,options.bf_channel);
            
            if options.ask_me
                im_combination = [padarray(im(:,:,:,options.bf_channel),[size(im_rotate,1) - size(im,1), size(im_rotate,2) - size(im,2)],'post'), im_rotate(:,:,:,options.bf_channel)];

                answer = cw.plot.confirm_results(im_combination(:,:,:,options.bf_channel),'Rotation results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end
            
            im = im_rotate;
        end
        
        % Save movie to .mat for quicker access times
        
        disp('Saving input movie to .mat...')
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
        
        if strcmp(options.processing_options.wseg_mode,'otsu')
            well_segmentation_results_struct = cw.process.segment_wells_otsu( im(:,:,:,options.bf_channel), options );
        elseif strcmp(options.processing_options.wseg_mode,'edge')
            well_segmentation_results_struct = cw.process.segment_wells_edge( im(:,:,:,options.bf_channel), options );
        else
            error('Unrecognized well segmentation mode.')
        end
                
        if options.ask_me
            im_combination = [im(:,:,:,options.bf_channel), well_segmentation_results_struct.im_seg_final];

            answer = cw.plot.confirm_results(im_combination(:,:,:),'Well segmentation results.');

            if ~strcmp(answer,'Yes')
                return
            end
        end        
        
        disp('Saving well segmentation results...')
        
        save([full_filename '__analysis_results/well_segmentation.mat'],'well_segmentation_results_struct')
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Well-by-well detection of signal presence
    
    if propts.start_at <= 3 || ~exist([full_filename '__analysis_results/noise_detection.mat'],'file') 
        
        disp('Detecting wells and frames with signal in each channel...')

        signal_detection_results_struct = cw.process.detect_noise( well_tracking_results_struct, propts );
        
        disp('Saving signal detection results....')
        
        save([full_filename '__analysis_results/noise_detection.mat'],'signal_detection_results_struct')
    else
        disp('Loading previous signal detection data...')
        
        signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell segmentation
        
    if propts.start_at <= 4 || ~exist([full_filename '__analysis_results/cell_segmentation.mat'],'file')
        
        disp('Segmenting cells in each well...')

        cell_segmentation_results_struct = cw.process.segment_cells( well_tracking_results_struct, ...
            signal_detection_results_struct, propts );
        
        disp('Cleaning up segmentation....')
        
        cell_segmentation_results_struct = cw.process.clean_up_cell_segmentation( cell_segmentation_results_struct, [] );
       
        disp('Saving cell segmentation results....')
        
        save([full_filename '__analysis_results/cell_segmentation.mat'],'cell_segmentation_results_struct')

    else
        disp('Loading previous cell segmentation data...')
        
        cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell tracking
    
    if propts.start_at <= 5 || ~exist([full_filename '__analysis_results/cell_tracking.mat'],'file')
        
        disp('Tracking cells in each well...')
        
        cell_tracking_results_struct = cw.process.track_cells( cell_segmentation_results_struct, propts );
        
        disp('Saving cell tracking results....')
        
        save([full_filename '__analysis_results/cell_tracking.mat'],'cell_tracking_results_struct')

    else
%         disp('Loading previous cell tracking...')
        
%         cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    end
end




