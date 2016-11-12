function cell_segmentation_results_struct_out = segment_cells_nuclei_assignment( well_tracking_results_struct, cell_segmentation_results_struct, options )
    propts = options.processing_options;
    
    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    
    cell_segmentation_results_struct_out = cell_segmentation_results_struct;
    
    for well_idx = 1:num_wells
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        cell_masks = zeros(size(cell_segmentation_results_struct.cell_masks{well_idx}));
        
        for frame_idx = 1:num_frames
            cur_slice = squeeze(cur_well_im(:,:,frame_idx,:));
            cur_nuclei_mask = cell_segmentation_results_struct.cell_masks{well_idx}(:,:,frame_idx,options.nuclei_channel);
                        
            %%% assign the nuclei objects to the cell channel they belong
            %%% to using Mann-Whitney statistical test

            nuclei_objects = cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,options.nuclei_channel};
            
            for nuclei_object_idx = 1:numel(nuclei_objects)
                % go through all objects in the nuclei channel and create a
                % new mask for every cell channel
                
                start_i = ceil(nuclei_objects(nuclei_object_idx).BoundingBox(2));
                end_i = ceil(start_i + nuclei_objects(nuclei_object_idx).BoundingBox(4)-1);
                
                start_j = ceil(nuclei_objects(nuclei_object_idx).BoundingBox(1));
                end_j = ceil(start_j + nuclei_objects(nuclei_object_idx).BoundingBox(3)-1);
                
                nucleus_image = nuclei_objects(nuclei_object_idx).FilledImage;
                
                % the background we will compare the signal to in each cell channel is everything that's not in the entire nucleus mask
                
                % which channel has a greater probability of being above
                % the background noise for this nucleus object?
                
                cell_channel_p_values = [];
                
                for cell_channel_idx = options.cell_channels
                    channel_slice = cur_slice(:,:,cell_channel_idx);
                    
                    % get the bitmap image of this cell channel after overlaying the nucleus mask
                    
                    cur_cell_channel_object_pixels = channel_slice(start_i:end_i,start_j:end_j);
                    cur_cell_channel_object_pixels(nucleus_image == 0) = 0;
                    
                    % get the background which is everything that wasn't determined to be a nucleus
                    
                    channel_slice_background = channel_slice;
                    channel_slice_background(cur_nuclei_mask(:) ~= 0) = 0;
                    channel_slice_background = reshape(channel_slice_background,size(cur_slice(:,:,cell_channel_idx)));

                    % calculate Mann-Whitney of signal compared to background
                    
                    p_value = ranksum(channel_slice_background(channel_slice_background ~= 0),cur_cell_channel_object_pixels(cur_cell_channel_object_pixels ~= 0));
                    
                    cell_channel_p_values = [cell_channel_p_values, p_value];
                end
                
                % the assignment is the channel with the lowest p-value
                % put this nucleus mask as the cell object mask
                
                CV = std(cell_channel_p_values)/mean(cell_channel_p_values);
                
                % this is to account for when the background is really bad in both channels and we might as well discard
                % CV appears to hover around 1.4 when P-values are highly distinuishable and there is real signal
                
                if CV > 1.2 
                    [~,cell_channel_assignment] = min(cell_channel_p_values);
                    cell_channel_assignment = cell_channel_assignment + options.cell_channels(1) - 1;

                    assigned_cell_mask = cell_masks(:,:,frame_idx,cell_channel_assignment);
                    assigned_cell_mask(start_i:end_i,start_j:end_j) = nucleus_image;
                    cell_masks(:,:,frame_idx,cell_channel_assignment) = assigned_cell_mask;
                end
            end
            
            for cell_channel_idx = options.cell_channels
            
                % now detect channel-assigned cell objects in cell channels

                im_thresh_final = bwlabeln(cell_masks(:,:,frame_idx,cell_channel_idx),4);
                objects_final = regionprops(im_thresh_final,'ConvexHull','Centroid','FilledImage','BoundingBox','Extrema','Area');

                cell_segmentation_results_struct_out.detected_cell_props{well_idx,frame_idx,cell_channel_idx} = objects_final';

                %%%%%%%%%%%%%%%
                % add a border to all the objects for visualization

                im_thresh_temp = im_thresh_final;
                im_thresh_temp(im_thresh_temp(:) ~= 0) = im_thresh_temp(im_thresh_temp(:) ~= 0) + 1;

                im_thresh_final = reshape(im_thresh_temp,size(im_thresh_final));

                for object_idx = 1:numel(objects_final)

                    start_i = min(ceil(objects_final(object_idx).Extrema(:,2)));
                    start_j = min(ceil(objects_final(object_idx).Extrema(:,1)));

                    imloc_i = [start_i:start_i + size(objects_final(object_idx).FilledImage,1)-1];
                    imloc_j = [start_j:start_j + size(objects_final(object_idx).FilledImage,2)-1];

                    border_image = bwperim(objects_final(object_idx).FilledImage);

                    ti = im_thresh_final(imloc_i,imloc_j);
                    ti2 = ti;
                    ti2(border_image == 1) = 1;
                    ti = reshape(ti2, size(ti));

                    im_thresh_final(imloc_i,imloc_j) = ti;
                end
                
                cell_masks(:,:,frame_idx,cell_channel_idx) = im_thresh_final;
            end
            
            cell_masks(:,:,frame_idx,options.nuclei_channel) = cur_nuclei_mask;
            
            if 0
                %%%% DEBUG PLOTTING
                
                nuclei_pixels = cur_slice(:,:,options.nuclei_channel);
                nuclei_pixels_masked = nuclei_pixels;
                nuclei_pixels_masked(cur_nuclei_mask == 0) = 0;
                nuclei_pixels_masked = reshape(nuclei_pixels_masked,size(nuclei_pixels));

                channel_1_mask = cell_masks(:,:,frame_idx,options.cell_channels(1));
                channel_1_pixels = cur_slice(:,:,options.cell_channels(1));
                channel_1_pixels_masked = channel_1_pixels;
                channel_1_pixels_masked(channel_1_mask == 0) = 0;
                channel_1_pixels_masked = reshape(channel_1_pixels_masked,size(channel_1_mask));

                channel_2_mask = cell_masks(:,:,frame_idx,options.cell_channels(2));
                channel_2_pixels = cur_slice(:,:,options.cell_channels(2));
                channel_2_pixels_masked = channel_2_pixels;
                channel_2_pixels_masked(channel_2_mask == 0) = 0;
                channel_2_pixels_masked = reshape(channel_2_pixels_masked,size(channel_2_mask));
                
                figure(123089)
                clf

                    subtightplot(3,3,1)
                        hold all

                        imagesc(nuclei_pixels)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                    subtightplot(3,3,4)
                        hold all

                        imagesc(nuclei_pixels_masked)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])
                        
                    subtightplot(3,3,7)
                        hold all

                        imagesc(cur_nuclei_mask)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                    %%%

                    subtightplot(3,3,2)
                        hold all

                        imagesc(channel_1_pixels)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                    subtightplot(3,3,5)
                        hold all

                        imagesc(channel_1_pixels_masked)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])
                        
                    subtightplot(3,3,8)
                        hold all

                        imagesc(channel_1_mask)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                    %%%

                    subtightplot(3,3,3)
                        hold all

                        imagesc(channel_2_pixels)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                    subtightplot(3,3,6)
                        hold all

                        imagesc(channel_2_pixels_masked)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])
                        
                    subtightplot(3,3,9)
                        hold all

                        imagesc(channel_2_mask)

                        axis image
                        set(gca,'Ydir','Reverse')

                        set(gca,'Color','white')
                        set(gca,'XTick',[])
                        set(gca,'YTick',[])

                suptitle(['Well: ' num2str(well_idx) ' Frame: ' num2str(frame_idx)])
                colormap gray
                %%%%
                
                pause
            end
            
        end
        
        cell_segmentation_results_struct_out.cell_masks{well_idx} = cell_masks;
    end
end

