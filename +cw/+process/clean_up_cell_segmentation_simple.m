function [cell_segmentation_results_corrected_struct,validation_images] = clean_up_cell_segmentation( cell_segmentation_results_struct, cell_quality_metrics )
    cell_quality_metrics.min_area = -inf; %pixels
    cell_quality_metrics.max_area = inf; %pixels
    
    object_mask_overlap_threshold = 10; %pixels
    
    % pull variables needed
    
    cell_segmentation_results_corrected_struct = cell_segmentation_results_struct;
    
    masks_water = cell_segmentation_results_struct.cell_masks_water;
    
    objects_nowater = cell_segmentation_results_struct.detected_cell_props_nowater;
    objects_water = cell_segmentation_results_struct.detected_cell_props_water;
    
    num_wells = size(objects_nowater,1);
    num_frames = size(objects_nowater,2);
    num_channels = size(objects_nowater,3);
    
    masks_final = masks_water;
    objects_final = cell(size(objects_nowater));
    
    % iterate through frames discarding bad objects as well as combining
    % watershedded and non-watershedded object maps when appropriate
    
    multiWaitbar('CloseAll');
    multiWaitbar('Cleaning up cell segmentation...',0);
    
    for well_idx = 1:num_wells
        for frame_idx = 1:num_frames

            %%% Remove objects that do not meet cell criteria

            for channel_idx = 1:num_channels
                objects_nowater{well_idx,frame_idx,channel_idx} = screen_objects(objects_nowater{well_idx,frame_idx,channel_idx}, cell_quality_metrics);
                objects_water{well_idx,frame_idx,channel_idx} = screen_objects(objects_water{well_idx,frame_idx,channel_idx}, cell_quality_metrics);
            end

            %%% Now make sure objects in the non-watershedded maps exist in
            %%% some form in the watershedded maps
            
            % in regular region properties, background is 0. all other objects are
            % numbered

            % in watershedding, background is 1, object perimeters are 0, objects
            % are numbered > 1
            
            for channel_idx = 1:size(objects_nowater,3)
%                 cur_mask_nowater = masks_nowater{well_idx}(:,:,frame_idx,channel_idx);
                cur_mask_water = masks_water{well_idx}(:,:,frame_idx,channel_idx);
                cur_corrected_mask = cur_mask_water;
                
                cur_nowater_objects = objects_nowater{well_idx,frame_idx,channel_idx};
                cur_water_objects = objects_water{well_idx,frame_idx,channel_idx};
                cur_objects_final = cur_water_objects;
                
                for obj_nowater_idx = 1:numel(cur_nowater_objects)
                    % Can we find an object in the watershedded image that
                    % overlaps with the non-watershedded object?
                    
                    % get pixel list of current nonwater object
                    
                    nowater_pixels = cur_nowater_objects(obj_nowater_idx).PixelIdxList;
                    
                    dont_add_nowater_obj = 0;
                    
                    for obj_water_idx = 1:numel(cur_water_objects)
                        % put all checked watershedded objs into list
                        
                        cur_objects_final = [cur_objects_final, cur_water_objects(obj_water_idx)];
                        
                        % get pixel list of current water object
                        
                        water_pixels = cur_water_objects(obj_water_idx).PixelIdxList;
                        
                        % does the current nonwater object overlap with
                        % this water object by at least
                        % object_mask_overlap_threshold pixels?
                        
                        object_overlap_area = sum(ismember(nowater_pixels,water_pixels));
                        
                        if object_overlap_area > object_mask_overlap_threshold
                            % these are two apparently overlapping objects.
                            % this excluded nonwater from being added
                            
                            dont_add_nowater_obj = 1;
                            
                            break 
                        end
                    end
                    
                    if ~dont_add_nowater_obj
                        % This nonwater object did not correspond to a
                        % watershedded obj so we should added it to the
                        % final mask
                        
                        next_obj_id = max(cur_corrected_mask(:))+1;
                        
                        cur_corrected_mask(nowater_pixels) = next_obj_id;
                        
                        % add a perimeter to the image
                        
                        perim_im_size_matched = zeros(size(cur_corrected_mask));
                        
                        perim_im = bwperim(cur_nowater_objects(obj_nowater_idx).Image);
                        
                        bb = round(cur_nowater_objects(obj_nowater_idx).BoundingBox);
                        
                        perim_im_size_matched(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = perim_im;
                        
%                         figure(314)
%                         clf
%                             subplot(1,2,1)
%                                 imagesc(perim_im)
%                                 axis image
%                             
%                             subplot(1,2,2)
%                                 imagesc(perim_im_size_matched)
%                                 axis image
%                         
%                         return
    
                        cur_corrected_mask(logical(perim_im_size_matched)) = 0;
                        cur_objects_final = [cur_objects_final, cur_nowater_objects(obj_nowater_idx)];
                        
                    end
                end
                
                masks_final{well_idx}(:,:,frame_idx,channel_idx) = cur_corrected_mask;
                objects_final{well_idx,frame_idx,channel_idx} = cur_objects_final;
            end
        end
        
        multiWaitbar('Cleaning up cell segmentation...',well_idx / num_wells);
    end
    
    cell_segmentation_results_corrected_struct.cell_masks_final = masks_final;
    cell_segmentation_results_corrected_struct.detected_cell_props_final = objects_final;
   
    multiWaitbar('CloseAll');
end

function good_cells = screen_objects(objects, cell_quality_metrics)
    good_cells = [];

    for obj_idx = 1:numel(objects)
        if cell_quality_metrics.min_area > objects(obj_idx).Area
            continue
        end
        
        if cell_quality_metrics.max_area < objects(obj_idx).Area
            continue
        end
        
        good_cells = [good_cells, objects(obj_idx)];
    end
end