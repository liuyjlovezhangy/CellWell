function [well_segmentation_results_struct, im_seg_final] = segment_wells_otsu2(im, options)

    im_bf = im(:,:,:,options.bf_channel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    well_well_dist_width = options.well_width + options.well_spacing_width;
    well_array_width = options.well_width * options.well_counts(2) + options.well_spacing_width * (options.well_counts(2)-1);
    
    well_well_dist_height = options.well_height + options.well_spacing_height;
    well_array_height = options.well_height * options.well_counts(1) + options.well_spacing_height * (options.well_counts(1)-1);

    min_well_area = 0.8 * options.well_width * options.well_height;
    max_well_area = 1.2 * options.well_width * options.well_height;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% init
    
    num_frames = size(im,3);
    num_channels = size(im,4);
    
    extrema_extrema_x = repmat([inf,1],[num_frames,1]);
    extrema_extrema_y = repmat([inf,1],[num_frames,1]);
    
    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run segmentation

    im_noise = zeros(size(im_bf));
    im_seg = zeros(size(im_bf));
    
    multiWaitbar('CloseAll');
    multiWaitbar('Segmenting...',0);
    
    for frame_idx = 1:num_frames
        
        cur_frame = im_bf(:,:,frame_idx);
        
        im_noise(:,:,frame_idx) = bpass(cur_frame,1,30);

        in = im_noise(:,:,frame_idx);

        seclose = strel('rect',[5 5]);
        seopen = strel('rect',[5 5]);
        
        im_seg(:,:,frame_idx) = imfill(imopen(imclose(imclearborder(in==0),seclose),seopen),'holes');

        if options.processing_options.wseg_debug
            figure(14332)
            clf
            
                subtightplot(1,3,1)
                    hold all
                    
                    imagesc(cur_frame)
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(1,3,2)
                    hold all
                    
                    imagesc(im_noise(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(1,3,3)
                    hold all
                    
                    imagesc(im_seg(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                pause
        end
        
        multiWaitbar('Segmenting...',frame_idx/num_frames);
    end
    
    im_mask = im_seg;
    
    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% infer the surrounding bounding box for all wells put together
    
    multiWaitbar('CloseAll');
    multiWaitbar('Detecting objects...',0);
    
    objects = cell(1,num_frames);
    
    %%% object detection to find extrema of well spacing in this frame
    for frame_idx = 1:num_frames
        im_seg(:,:,frame_idx) = im_bf(:,:,frame_idx) .* im_mask(:,:,frame_idx);
        
        L = bwlabeln(im_seg(:,:,frame_idx));
        objects{frame_idx} = regionprops(L,im_seg(:,:,frame_idx),'Extrema','Area','Centroid');
        
        for obj_idx = 1:numel(objects{frame_idx})
            if objects{frame_idx}(obj_idx).Area > min_well_area && objects{frame_idx}(obj_idx).Area < max_well_area
                
                extrema_extrema_x(frame_idx,1) = min([extrema_extrema_x(frame_idx,1), min(objects{frame_idx}(obj_idx).Extrema(:,1))]);
                extrema_extrema_x(frame_idx,2) = max([extrema_extrema_x(frame_idx,2), max(objects{frame_idx}(obj_idx).Extrema(:,1))]);
                
                extrema_extrema_y(frame_idx,1) = min([extrema_extrema_y(frame_idx,1), min(objects{frame_idx}(obj_idx).Extrema(:,2))]);
                extrema_extrema_y(frame_idx,2) = max([extrema_extrema_y(frame_idx,2), max(objects{frame_idx}(obj_idx).Extrema(:,2))]);
            end

        end
        
        multiWaitbar('Detecting objects...',frame_idx/num_frames);
    end
    
    
    
    multiWaitbar('CloseAll');
    
    mask_final = zeros(size(im_bf));
    im_seg_final = zeros(size(im));
    
    %%% repair extrema matrices for frames that didn't segment well by
    %%% borrowing extrema from neighboring frames
    
    % sanity check
    
    extrema_extrema_x(extrema_extrema_x(:,1) + well_array_width > size(im_bf,2),:) = NaN;
    extrema_extrema_y(extrema_extrema_y(:,1) + well_array_height > size(im_bf,1),:) = NaN;
    
    mask_width = extrema_extrema_x(:,2) - extrema_extrema_x(:,1);
    mask_height = extrema_extrema_y(:,2) - extrema_extrema_y(:,1);
    
    extrema_extrema_x(mask_width < 0.9*well_array_width | mask_width > 1.1*well_array_width,:) = NaN;
    extrema_extrema_y(mask_height < 0.9*well_array_height | mask_height > 1.1*well_array_height,:) = NaN;
    
    % outlier detection
    
    alpha = 0.25;
    
    extrema_extrema_x = [deleteoutliers(extrema_extrema_x(:,1),alpha,1), deleteoutliers(extrema_extrema_x(:,2),alpha,1)];
    extrema_extrema_y = [deleteoutliers(extrema_extrema_y(:,1),alpha,1), deleteoutliers(extrema_extrema_y(:,2),alpha,1)];
    
    % sweep extrema list and set NaNs to nearest neighbor extrema
    
    for frame_idx = 2:num_frames
        if isnan(sum(extrema_extrema_x(frame_idx,:))) || isnan(sum(extrema_extrema_y(frame_idx,:)))
            extrema_extrema_x(frame_idx,:) = extrema_extrema_x(frame_idx-1,:);
            extrema_extrema_y(frame_idx,:) = extrema_extrema_y(frame_idx-1,:);
        end
    end
    
    for frame_idx = num_frames-1:-1:1
        if isnan(sum(extrema_extrema_x(frame_idx,:))) || isnan(sum(extrema_extrema_y(frame_idx,:)))
            extrema_extrema_x(frame_idx,:) = extrema_extrema_x(frame_idx+1,:);
            extrema_extrema_y(frame_idx,:) = extrema_extrema_y(frame_idx+1,:);
        end
    end
    
    % make frame-by-frame masks
    
    for frame_idx = 1:num_frames
        
        for row_idx = 1:options.well_counts(1)
            for col_idx = 1:options.well_counts(2)
                start_i = round(extrema_extrema_y(frame_idx,1)) + (row_idx-1) * well_well_dist_height;
                end_i = start_i + options.well_height;
                start_j = round(extrema_extrema_x(frame_idx,1)) + (col_idx-1) * well_well_dist_width + options.processing_options.wseg_otsu2_left_offset;
                end_j = start_j + options.well_width;

                mask_final(start_i:end_i, start_j:end_j, frame_idx) = 1;
            end 
        end
        
        for channel_idx = 1:num_channels
            im_seg_final(:,:,frame_idx,channel_idx) = im(:,:,frame_idx,channel_idx) .* mask_final(:,:,frame_idx);
        end
    end
    
    % detect actual wells based on our inferred mask
    
    for frame_idx = 1:num_frames
        L = bwlabeln(mask_final(:,:,frame_idx));
        objects = regionprops(L,'Centroid','Extrema');
        
        if numel(objects) ~= prod(options.well_counts)
            error('The segmentation did not work properly...')
        end
        
        well_segmentation_results_struct.frame(frame_idx).all_objects = objects;
        well_segmentation_results_struct.frame(frame_idx).good_objects = objects;
    end
    
%     well_segmentation_results_struct.im_seg_final = im_seg_final;
end
