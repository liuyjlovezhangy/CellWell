function [cell_segmentation_results_struct,validation_images] = segment_cells_simple( well_tracking_results_struct, signal_detection_results_struct, options )

    error('segment_cells_simple is broken. fix me')

    propts = options.processing_options;

    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    watershed_inputs = cell(1,num_wells);
    
    cell_masks_nowater = cell(1,num_wells);
    cell_masks_water = cell(1,num_wells);
    
    detected_cell_props_nowater = cell(num_wells,num_frames,num_channels-1);
    detected_cell_props_water = cell(num_wells,num_frames,num_channels-1);
    
    all_threshold_levels = cell(1,num_wells);
    all_thresh_xvals = cell(1,num_wells);
    all_thresh_yvals = cell(1,num_wells);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Performing cell segmentation...',0);
    
    validation_images = cell(1,(num_channels-1) * num_wells);
    
    for well_idx = 1:num_wells
        
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        cur_well_im_thresh_nowater = zeros(size(cur_well_im,1),size(cur_well_im,2),size(cur_well_im,3),num_channels-1);
        cur_well_im_thresh_water = zeros(size(cur_well_im,1),size(cur_well_im,2),size(cur_well_im,3),num_channels-1);
        cur_watershed_inputs = zeros(size(cur_well_im,1),size(cur_well_im,2),size(cur_well_im,3),num_channels-1);

        cur_threshold_levels = cell(size(cur_well_im,3),num_channels-1);

        cur_thresh_xvals = cell(size(cur_well_im,3),num_channels-1);
        cur_thresh_yvals = cell(size(cur_well_im,3),num_channels-1);

        multiWaitbar('Current frame...',0);

        for frame_idx = 1:num_frames

            for channel_idx = 1:num_channels-1
                L_nowater = zeros(size(cur_well_im,1),size(cur_well_im,2));
                L_water = zeros(size(cur_well_im,1),size(cur_well_im,2));
                D = zeros(size(cur_well_im,1),size(cur_well_im,2));
                
                if signal_detection_results_struct.is_noise_matrix(well_idx,channel_idx,frame_idx)
                    cur_cell_mask = zeros(size(cur_well_im,1),size(cur_well_im,2));
                else
                    if frame_idx > 1 
                        levels_prev = cur_threshold_levels{frame_idx-1,channel_idx};
                    else
                        levels_prev = [];
                    end
                    
                    [cur_threshold_levels{frame_idx,channel_idx}.level, cur_threshold_levels{frame_idx,channel_idx}.level_low, cur_threshold_levels{frame_idx,channel_idx}.level_high, ...
                        cur_thresh_xvals{frame_idx,channel_idx}, cur_thresh_yvals{frame_idx,channel_idx}] = ...
                        cw.process.adaptiveThresholdFinder(cur_well_im(:,:,frame_idx,channel_idx), propts.cseg_adaptive_thresh_scale,...
                        propts.cseg_tracking_params, levels_prev);

                    cur_cell_mask = threshold3D(cur_well_im(:,:,frame_idx,channel_idx), cur_threshold_levels{frame_idx,channel_idx}.level);

                    %%% Image morphological operations to clean up

                    % link together objects

                    se = strel('disk',propts.cseg_close_radius(channel_idx));
                    cur_cell_mask = imclose(cur_cell_mask,se);
                    
                    % Delete very small connections between cells
                    
                    se = strel('disk',propts.cseg_open_radius(channel_idx));
                    cur_cell_mask = imopen(cur_cell_mask,se);

                    % Clear very small objects
                    cur_cell_mask = bwareaopen(cur_cell_mask, propts.cseg_min_cell_area);

                    im_seg = cur_well_im(:,:,frame_idx,channel_idx) .* cur_cell_mask;

                    if propts.cseg_watershedding(channel_idx)
                        %%% watershed to split up close cells

%                         im_seg_water = adapthisteq(im_seg);
                        im_seg_water = im_seg;

                        if strcmp(propts.cseg_water_method,'dist')
                        
                            D = bwdist(~cur_cell_mask);
                            D = -D;
                            
                        elseif strcmp(propts.cseg_water_method,'intensity')
                            cell_mask_em = imextendedmax(im_seg_water, 0.9);

                            % D is the input to watershed after the initial
                            % segmented image is cleaned up

                            D = im_seg_water;

                            h = fspecial('gaussian',[5 5], 10);
                            D = imfilter(D,h, 'replicate');

                            D = imcomplement(D);

                            D = imimposemin(D, ~cur_cell_mask | cell_mask_em);
                        else
                            error('Given watershed method does not exist')
                        end
                        
                        D(~cur_cell_mask) = -inf;
                        
                        
                        
                        L_water = watershed(D);
                        
                        

                        % Object detection for watershed mask

                        water_objects = regionprops(L_water,im_seg_water,'all');

                        % The largest area in this image will be the
                        % background unless something is very screwed up or
                        % you have loaded a carpet of cells
                        
                        water_areas = [water_objects.Area];
                        [~,background_idx] = max(water_areas);
                        water_objects(background_idx) = [];

                        detected_cell_props_water{well_idx,frame_idx,channel_idx} = water_objects;    
                    end
                    
                    L_nowater = bwlabeln(cur_cell_mask,4);
                    
                    % Object detection for non-watershed mask

                    nowater_objects = regionprops(L_nowater,im_seg,'all');

                    detected_cell_props_nowater{well_idx,frame_idx,channel_idx} = nowater_objects;
                end

                cur_well_im_thresh_nowater(:,:,frame_idx,channel_idx) = L_nowater;
                cur_well_im_thresh_water(:,:,frame_idx,channel_idx) = L_water;
                cur_watershed_inputs(:,:,frame_idx,channel_idx) = D;
                
                
            end

            multiWaitbar('Current frame...',frame_idx / num_frames);
            
            
        end

        cell_masks_nowater{well_idx} = cur_well_im_thresh_nowater;
        cell_masks_water{well_idx} = cur_well_im_thresh_water;
        watershed_inputs{well_idx} = cur_watershed_inputs;
        
        all_threshold_levels{well_idx} = cur_threshold_levels;
        all_thresh_xvals{well_idx} = cur_thresh_xvals;
        all_thresh_yvals{well_idx} = cur_thresh_yvals;

        multiWaitbar('Performing cell segmentation...',well_idx / num_wells);
    end
    
    if options.ask_me
        % build a validation image array
        
        for well_idx = 1:num_wells
            for channel_idx = 1:num_channels-1
                cur_well_im =  mat2gray(well_tracking_results_struct.wells(well_idx).im_well(:,:,:,channel_idx));
                cur_well_im_thresh_nowater = cell_masks_nowater{well_idx}(:,:,:,channel_idx);
                cur_watershed_inputs = watershed_inputs{well_idx}(:,:,:,channel_idx);
                cur_well_im_thresh_water = cell_masks_water{well_idx}(:,:,:,channel_idx);
                
                valim = [cur_well_im, cur_well_im_thresh_nowater, cur_watershed_inputs, cur_well_im_thresh_water];
                
                validation_images{sub2ind([num_channels-1,num_wells],channel_idx,well_idx)} = valim;
            end
        end
        
    end
    
    
    cell_segmentation_results_struct.cell_masks_nowater = cell_masks_nowater;
    cell_segmentation_results_struct.cell_masks_water = cell_masks_water;
    cell_segmentation_results_struct.watershed_inputs = watershed_inputs;
    cell_segmentation_results_struct.threshold_levels = all_threshold_levels;
    cell_segmentation_results_struct.thresh_xvals = all_thresh_xvals;
    cell_segmentation_results_struct.thresh_yvals = all_thresh_yvals;
    cell_segmentation_results_struct.detected_cell_props_nowater = detected_cell_props_nowater;
    cell_segmentation_results_struct.detected_cell_props_water = detected_cell_props_water;
    cell_segmentation_results_struct.cell_segmentation_opts = propts;
    
    multiWaitbar('CloseAll');
    drawnow
end
