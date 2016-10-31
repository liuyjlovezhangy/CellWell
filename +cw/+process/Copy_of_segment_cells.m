function cell_segmentation_results_struct = segment_cells( well_tracking_results_struct, signal_detection_results_struct, options )

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
                        cw.process.adaptiveThresholdFinder(cur_well_im(:,:,frame_idx,channel_idx), options.cseg_adaptive_thresh_scale,...
                        options.cseg_tracking_params, levels_prev);

                    cur_cell_mask = threshold3D(cur_well_im(:,:,frame_idx,channel_idx), cur_threshold_levels{frame_idx,channel_idx}.level);

                    %%% Image morphological operations to clean up

                    % link together objects

                    se = strel('disk',options.cseg_close_radius(channel_idx));
                    cur_cell_mask = imclose(cur_cell_mask,se);

                    % Clear very small objects
                    cur_cell_mask = bwareaopen(cur_cell_mask, options.cseg_min_cell_area);

                    im_seg = cur_well_im(:,:,frame_idx,channel_idx) .* cur_cell_mask;

                    if options.cseg_watershedding(channel_idx)
                        % watershed to split up close cells
                        
%                         D = bwdist(~cur_cell_mask);
%                         D = -D;

                        im_seg_water = adapthisteq(im_seg);
                        cell_mask_perim = bwperim(cur_cell_mask);

                        cell_mask_em = imextendedmax(im_seg_water, 0.3);

                        % D is the input to watershed after the initial
                        % segmented image is cleaned up
                        
                        D = im_seg_water;
                        
%                         figure
%                             imagesc(D)
                        
                        h = fspecial('gaussian',[5 5], 10);
                        D = imfilter(D,h, 'replicate');
                        
%                         figure
%                             imagesc(D)
%                         error
                        
                        D = imcomplement(D);
                        
                        D = imimposemin(D, ~cur_cell_mask | cell_mask_em);
                        
                        
                        
                        D(~cur_cell_mask) = -inf;
                        
                        if 0
                            figure(123423)
                            clf
                                subtightplot(1,4,1)
                                    imshow(~cur_cell_mask)
                                    axis image
                                subtightplot(1,4,2)
                                    imshow(cell_mask_em)
                                    axis image
                                subtightplot(1,4,3)
                                    imshow(~cur_cell_mask | cell_mask_em)
                                    axis image
                                subtightplot(1,4,4)
                                    imshow(D)
                                    axis image
                                    
%                                     error
                        end


                        L_water = watershed(D);

                        if 0
                            figure(13412)
                            clf
                                subplot(2,3,1)
                                    imagesc(im_seg)
                                    axis image
                                subplot(2,3,2)
                                    imagesc(cur_cell_mask)
                                    axis image
                                subplot(2,3,3)
                                    A=imoverlay(im_seg, cell_mask_perim, [.3 1 .3]);
                                    imagesc(A)
                                    axis image
                                subplot(2,3,4)
                                    imagesc(D)
                                    axis image
                                subplot(2,3,5)
                                    imagesc(L)
                                    axis image
                                subplot(2,3,6)
                                    A=imoverlay(im_seg, ~L, [.3 1 .3]);
                                    imagesc(A)
                                    axis image
                            error
                        end
                        
                        if 0
                            figure(1034)
                            clf

                                subplot(2,3,1)
                                hold all

                                    imagesc(cur_well_im(:,:,frame_idx,channel_idx))

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Well image')

                                subplot(2,3,2)
                                hold all

                                    imagesc(cur_cell_mask)

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Well mask')

                                subplot(2,3,3)
                                hold all

                                    imagesc(im_seg)

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Well image masked')

                                subplot(2,3,4)
                                hold all

                                    imagesc(L_regular)

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Well object regular detection')

                                subplot(2,3,5)
                                hold all

                                    imagesc(D)

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Input to watershed')

                                subplot(2,3,6)
                                hold all

                                    imagesc(L_water)

                                    axis image
                                    axis tight
                                    set(gca,'Ydir','Reverse')
                                    axis off
                                    colormap gray

                                    title('Well object watershedding')

                            suptitle('Watershedding experiment')

                            set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
                            set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',3)
                            set(gcf, 'color', 'white');

                            drawnow

                            error

                        end
                        
                        % Object detection for watershedded cells in this channel

                        objects = regionprops(L_water,im_seg_water,'all');

                        areas = [objects.Area];

                        [~,background_idx] = max(areas);
                        objects(background_idx) = [];

                        good_objects = screen_good_objects(objects);
                        detected_cell_props_water{well_idx,frame_idx,channel_idx} = good_objects;    
                    end
                    
                    L_nowater = bwlabeln(cur_cell_mask);
                    
                    % Object detection for watershedded cells in this channel

                    objects = regionprops(L_nowater,im_seg,'all');
                    good_objects = screen_good_objects(objects);
                    detected_cell_props_nowater{well_idx,frame_idx,channel_idx} = good_objects;
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
    
    cell_segmentation_results_struct.cell_masks_nowater = cell_masks_nowater;
    cell_segmentation_results_struct.cell_masks_water = cell_masks_water;
    cell_segmentation_results_struct.watershed_inputs = watershed_inputs;
    cell_segmentation_results_struct.threshold_levels = all_threshold_levels;
    cell_segmentation_results_struct.thresh_xvals = all_thresh_xvals;
    cell_segmentation_results_struct.thresh_yvals = all_thresh_yvals;
    cell_segmentation_results_struct.detected_cell_props_nowater = detected_cell_props_nowater;
    cell_segmentation_results_struct.detected_cell_props_water = detected_cell_props_water;
    cell_segmentation_results_struct.cell_segmentation_opts = options;
    
    multiWaitbar('CloseAll');
    drawnow
end

function good_objects = screen_good_objects(objects)
    good_objects = [];
    
    for obj_idx = 1:numel(objects)
        % screen for bad cells

        bad_object_flag = 0;

        % if good, put it in the list

        if ~bad_object_flag
            good_objects = [good_objects, objects(obj_idx)];
        end
    end
end

