function well_segmentation_results_struct = segment_wells( im_bf, o )

    num_frames = size(im_bf,3);
    
    well_segmentation_results_struct = [];

    multiWaitbar('CloseAll');
    multiWaitbar('Performing well segmentation...',0);

    for frame_idx = 1:num_frames

        im_bf_slice = im_bf(:,:,frame_idx);

        %%% Find otsu threshold levels and binarize

%         [IDX,~,well_segmentation_results_struct.otsu_idcs] = otsu(im_bf_slice,o.wseg_num_otsu_levels);
        IDX = otsu(im_bf_slice,o.wseg_num_otsu_levels);

        im_mask_thresh = IDX == o.wseg_num_otsu_levels;

        %%% Image morphological operations to clean up

        % remove objects on border

        im_mask_noborder = imclearborder(im_mask_thresh);

        % link together objects

        se = strel('disk',o.wseg_close_size);
        im_mask_connect = imclose(im_mask_noborder,se);

        % delete border objects again

        im_mask_noborder2 = imclearborder(im_mask_connect);

        % Clear very small objects
        im_mask_cleared = bwareaopen(im_mask_noborder2, o.wseg_min_size);

        im_mask_final = im_mask_cleared;
        im_mask_final_corrected = im_mask_final;

        im_seg = im_bf_slice .* im_mask_final;

        %%% Object detection

        L = bwlabeln(im_seg);
        objects = regionprops(L,im_seg,'all');

        good_objects = [];

        %%% Cull bad objects

        for obj_idx = 1:numel(objects)
            if cw.process.well_criterion(objects(obj_idx))
                good_objects = [good_objects, objects(obj_idx)];

                extrema = objects(obj_idx).Extrema;

                xmin = floor(min(extrema(:,1)) - o.wseg_extra_border_x);
                xmax = ceil(max(extrema(:,1)) + o.wseg_extra_border_x);
                ymin = floor(min(extrema(:,2)) - o.wseg_extra_border_y);
                ymax = ceil(max(extrema(:,2)) + o.wseg_extra_border_y);

                im_mask_final_corrected(ymin:ymax,xmin:xmax) = 1;
            end
        end

        well_segmentation_results_struct.frame(frame_idx).all_objects = objects;
        well_segmentation_results_struct.frame(frame_idx).good_objects = good_objects;
        well_segmentation_results_struct.frame(frame_idx).mask = im_mask_final;
        well_segmentation_results_struct.frame(frame_idx).im_mask_thresh = im_mask_thresh;
        well_segmentation_results_struct.frame(frame_idx).im_mask_final_corrected = im_mask_final_corrected;
        
        multiWaitbar('Performing well segmentation...', frame_idx / num_frames);
    end
    
    well_segmentation_results_struct.extra_border_x = o.wseg_extra_border_x;
    well_segmentation_results_struct.extra_border_y = o.wseg_extra_border_y;
    
    multiWaitbar('CloseAll');
    drawnow
end

