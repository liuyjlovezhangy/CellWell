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
        
        if options.last_frame ~= -1
            im = im(:,:,options.first_frame:options.last_frame,:);
        else
            im = im(:,:,options.first_frame:end,:);
        end
        
        % Perform registration and rotation if desired
        
        if propts.register
            im_reg = cw.process.register(im,options.bf_channel);
            
            if options.ask_me

                im_combination = [im(:,:,:,options.bf_channel), im_reg(:,:,:,options.bf_channel)];

                answer = cw.plot.confirm_results(im_combination,'Registration results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end
            
            im = im_reg;
        end
        
        if propts.rotate
            im_rotate = cw.process.rotate(im,options.bf_channel);
            
            if options.ask_me
                if any([size(im_rotate,1) - size(im,1), size(im_rotate,2) - size(im,2)] < 0)
                    im_combination = [im(:,:,:,options.bf_channel), ...
                        padarray(im_rotate(:,:,:,options.bf_channel),[size(im,1) - size(im_rotate,1), size(im,2) - size(im_rotate,2)],'post')];
                else
                    im_combination = [padarray(im(:,:,:,options.bf_channel),[size(im_rotate,1) - size(im,1), size(im_rotate,2) - size(im,2)],'post'), im_rotate(:,:,:,options.bf_channel)];
                end

                answer = cw.plot.confirm_results(im_combination,'Rotation results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end
            
            im = im_rotate;
        end
        
        % Save movie to .mat for quicker access times
        
        disp('Saving input movie to .mat...')
        save([full_filename '__analysis_results/input_movie.mat'],'im','-v7.3')
        
    else
        disp('Loading movie from mat...')
        
        im = importdata([full_filename '__analysis_results/input_movie.mat']);
    end
    
    im = mat4D_to_gray(im);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform segmentation on each BF frame
    
    if propts.start_at <= 1 || ~exist([full_filename '__analysis_results/well_segmentation.mat'],'file')
    
        disp('Performing well segmentation...')
        
        if strcmp(propts.wseg_mode,'otsu')
            error('dont use this well segmentation mode')
            [well_segmentation_results_struct,im_seg_final] = cw.process.segment_wells_otsu( im, options );
        elseif strcmp(propts.wseg_mode,'otsu2')
            error('dont use this well segmentation mode')
            [well_segmentation_results_struct,im_seg_final] = cw.process.segment_wells_otsu2( im, options );
        elseif strcmp(propts.wseg_mode,'bandpass')
            error('dont use this well segmentation mode')
            [well_segmentation_results_struct,im_seg_final] = cw.process.segment_wells_bandpass( im, options );
        elseif strcmp(propts.wseg_mode,'template')
            [well_segmentation_results_struct,im_seg_final] = cw.process.segment_wells_template( im, options );
        else
            error('Unrecognized well segmentation mode.')
        end

        fresh_wseg = 1;
        
        disp('Saving well segmentation results...')
        
        save([full_filename '__analysis_results/well_segmentation.mat'],'well_segmentation_results_struct')
    else
        disp('Loading previous well segmentation...')
        
        fresh_wseg = 0;
        
        well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Perform well tracking to link wells from frame to frame
    %%% This generates super large files so now we just compute it every
    %%% time because it is computationally very fast.
    
    disp('Performing well tracking...')
        
    well_tracking_results_struct = cw.process.track_wells( im, well_segmentation_results_struct, options );
    
    if options.ask_me && fresh_wseg

        im_bf_all = cell(1,numel(well_tracking_results_struct.wells));
        
        for well_idx = 1:numel(well_tracking_results_struct.wells)
            im_bf_all{well_idx} = well_tracking_results_struct.wells(well_idx).im_well(:,:,:,options.bf_channel);
        end
        
        answer = cw.plot.confirm_results(im_bf_all,'Well segmentation results.');

        if ~strcmp(answer,'Yes')
            return
        end
    end       
    
    if options.export_tifs && fresh_wseg
        disp('Writing well tifs...')
        
        delete([full_filename '__analysis_results/well_im_*'])
        
        for well_idx = 1:numel(well_tracking_results_struct.wells)
            well_idx
            
            im = well_tracking_results_struct.wells(well_idx).im_well;
            
            zwriteim([full_filename '__analysis_results/well_im_' num2str(well_idx)],im);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Well-by-well detection of signal presence
    
    if propts.start_at <= 3 || ~exist([full_filename '__analysis_results/noise_detection.mat'],'file') 
        
        disp('Detecting wells and frames with signal in each channel...')

        [signal_detection_results_struct,validation_images] = cw.process.detect_noise( well_tracking_results_struct, options );
        
%         if options.ask_me
% 
%             answer = cw.plot.confirm_results(validation_images,'Noise detection results.');
% 
%             if ~strcmp(answer,'Yes')
%                 return
%             end
%         end
        
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

        if strcmp(propts.cseg_mode,'simple')
        
            error('dont use this cell segmentation mode')
            
            [cell_segmentation_results_struct,validation_images] = cw.process.segment_cells_simple( well_tracking_results_struct, ...
                signal_detection_results_struct, options );

            if options.ask_me

                answer = cw.plot.confirm_results(validation_images,'Cell segmentation results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end

            disp('Cleaning up segmentation....')

            [cell_segmentation_results_struct,validation_images] = cw.process.clean_up_cell_segmentation_simple( cell_segmentation_results_struct, [] );

            if options.ask_me

                answer = cw.plot.confirm_results(validation_images,'Cell segmentation CLEANING results.');

                if ~strcmp(answer,'Yes')
                    return
                end
            end

        elseif strcmp(propts.cseg_mode,'circle')
            error('dont use this cell segmentation mode')
            cell_segmentation_results_struct = cw.process.segment_cells_circle( well_tracking_results_struct, ...
                signal_detection_results_struct, options );
        elseif strcmp(propts.cseg_mode,'circle_bandpass')
            error('dont use this cell segmentation mode')
            cell_segmentation_results_struct = cw.process.segment_cells_circle_bandpass( well_tracking_results_struct, ...
                signal_detection_results_struct, options );
        elseif strcmp(propts.cseg_mode,'radial')
            cell_segmentation_results_struct = cw.process.segment_cells_radial_sym( well_tracking_results_struct, ...
                signal_detection_results_struct, options );
        elseif strcmp(propts.cseg_mode,'point_source')
            error('dont use this cell segmentation mode')
            cell_segmentation_results_struct = cw.process.segment_cells_point_source( well_tracking_results_struct, ...
                signal_detection_results_struct, options );    
        elseif strcmp(propts.cseg_mode,'mser')
            error('dont use this cell segmentation mode')
            cell_segmentation_results_struct = cw.process.segment_cells_mser( well_tracking_results_struct, ...
                signal_detection_results_struct, options );  
        else
            error('Unknown cell segmentation mode.')
        end
        
        if options.nuclei_channel
            cell_segmentation_results_struct = cw.process.segment_cells_nuclei_assignment(well_tracking_results_struct, cell_segmentation_results_struct, options);
        end
        
        disp('Saving cell segmentation results....')
        
        save([full_filename '__analysis_results/cell_segmentation.mat'],'cell_segmentation_results_struct')
        fresh_cseg = 1;

    else
        disp('Loading previous cell segmentation data...')
        
        cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
        
        fresh_cseg = 0;
    end

    if options.export_tifs && fresh_cseg
        disp('Writing cell segmentation tifs...')
        
        delete([full_filename '__analysis_results/cell_seg_im_*'])
        
        cell_masks = cell_segmentation_results_struct.cell_masks;
        
        for well_idx = 1:numel(cell_masks)
            well_idx
            
            im_well = mat4D_to_gray(well_tracking_results_struct.wells(well_idx).im_well);
            cur_mask = mat4D_to_gray(cell_masks{well_idx});
%             cur_mask(:,:,:,options.bf_channel) = [];
            
            zwriteim([full_filename '__analysis_results/cell_seg_im_' num2str(well_idx)],[im_well,cur_mask]);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell tracking
    
    if propts.start_at <= 5 || ~exist([full_filename '__analysis_results/cell_tracking.mat'],'file')
        
        disp('Tracking cells in each well...')
        
        cell_tracking_results_struct = cw.process.track_cells( cell_segmentation_results_struct, options );
        
        disp('Saving cell tracking results....')
        
        save([full_filename '__analysis_results/cell_tracking.mat'],'cell_tracking_results_struct')

    else
        disp('Loading previous cell tracking...')
        
        cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    end
    
%     delete([full_filename '__analysis_results/cell_tracks_well_*.csv'])
%     
%     for cell_channel_idx = 1:numel(options.tracking_channels)
%         
%         cur_tracks = cell_tracking_results_struct.cell_tracks(:,options.tracking_channels(cell_channel_idx));
%         
%         for well_idx = 1:numel(cur_tracks)
%             cur_cell = cur_tracks{well_idx};
%             
%             M = cur_cell;
%             if ~isempty(M)
%                 M = cellfun(@transpose,M,'UniformOutput',false);
%             
%                 num_tracks = numel(M);
%                 
%                 % make headers
%                 
%                 header = {};
%                 
%                 for track_idx = 1:num_tracks
% %                     header = [header, ['track_' num2str(track_idx) '_X']];
% %                     header = [header, ['track_' num2str(track_idx) '_Y']];
%                     
%                     header = [header, ['track_' num2str(track_idx)]];
%                 end
%                 
%                 M = [header;M]
%                 
% %                 dlmwrite([full_filename '__analysis_results/cell_tracks_well_' num2str(well_idx) '_channel_' num2str(options.tracking_channels(cell_channel_idx)) '.csv'],header,',');
% %                 
% %                 dlmwrite([full_filename '__analysis_results/cell_tracks_well_' num2str(well_idx) '_channel_' num2str(options.tracking_channels(cell_channel_idx)) '.csv'],M,'-append','delimiter',',');
%                 
%                 dlmcell([full_filename '__analysis_results/cell_tracks_well_' num2str(well_idx) '_channel_' num2str(options.tracking_channels(cell_channel_idx)) '.csv'],M,'delimiter',',');
%             end
%         end
%     end
%     
%     return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Interaction detection
    
    if propts.start_at <= 6 || ~exist([full_filename '__analysis_results/cell_interactions.mat'],'file')
        
        disp('Detecting cell interactions in each well...')
        
        cell_interaction_results_struct = cw.process.detect_interactions( well_tracking_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, options );
        
        disp('Saving cell interaction detection results....')
        
        save([full_filename '__analysis_results/cell_interactions.mat'],'cell_interaction_results_struct')

    else
        disp('Loading previous cell tracking...')
        
        cell_interaction_results_struct = importdata([full_filename '__analysis_results/cell_interactions.mat']);
    end
    
    return
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell death detection
    
    if propts.start_at <= 7 || ~exist([full_filename '__analysis_results/cell_interactions.mat'],'file')
        cell_death_results_struct = cw.process.detect_death( well_tracking_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, options );
        
        save([full_filename '__analysis_results/cell_death.mat'],'cell_death_results_struct')
    else
        disp('Loading previous cell death detection...')
        
        cell_interaction_results_struct = importdata([full_filename '__analysis_results/cell_death.mat']);
    end
end




