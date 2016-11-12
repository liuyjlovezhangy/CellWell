function cell_segmentation_results_struct_out = segment_cells_circle_nuclei_cleanup( well_tracking_results_struct, cell_segmentation_results_struct, options )
    propts = options.processing_options;
    
    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    cell_segmentation_results_struct_out = cell_segmentation_results_struct;
    
    for well_idx = 1:num_wells
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
%         well_masks = cell_segmentation_results_struct.cell_masks{well_idx};
        
        for frame_idx = 1:num_frames
            cur_slice = squeeze(cur_well_im(:,:,frame_idx,:));
%             cur_nuclei_mask = well_masks(:,:,frame_idx,options.nuclei_channel);
                        
            %%% assign the nuclei objects to the cell channel they belong
            %%% to using Mann-Whitney statistical test

            nuclei_objects = cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,options.nuclei_channel};
            
            for nuclei_object_idx = 1:numel(nuclei_objects)
                % go through all objects in the nuclei channel and create a
                % new mask for every cell channel
                
                start_i = floor(nuclei_objects(nuclei_object_idx).BoundingBox(2));
                end_i = ceil(start_i + nuclei_objects(nuclei_object_idx).BoundingBox(4));
                
                start_j = floor(nuclei_objects(nuclei_object_idx).BoundingBox(1));
                end_j = ceil(start_j + nuclei_objects(nuclei_object_idx).BoundingBox(3));
                
                % which channel has a greater probability of being above
                % the background noise for this nucleus object?
                
                
            end
            
            error
            
%             for cell_channel_idx = options.cell_channels
%                 channel_slice_inside = cur_slice(:,:,cell_channel_idx);
%                 channel_slice_outside = channel_slice_inside;
%                 
%                 channel_slice_inside(cur_nuclei_mask(:) == 0) = 0;
%                 channel_slice_outside(cur_nuclei_mask(:) ~= 0) = 0;
% 
%                 channel_slice_inside = reshape(channel_slice_inside,size(cur_slice(:,:,cell_channel_idx)));
%                 channel_slice_outside = reshape(channel_slice_outside,size(cur_slice(:,:,cell_channel_idx)));
%                 
%                 p_value = ranksum(channel_slice_inside(channel_slice_inside ~= 0),channel_slice_outside(channel_slice_outside ~= 0))
%                 
%                 figure(1342)
%                 clf
%                 
%                     subtightplot(1,2,1)
%                         hold all
%                         
%                         imagesc(channel_slice_inside)
% 
%                         axis image
%                         set(gca,'Ydir','Reverse')
% 
%                         set(gca,'Color','white')
%                         set(gca,'XTick',[])
%                         set(gca,'YTick',[])
%                         
%                     subtightplot(1,2,2)
%                         hold all
%                         
%                         imagesc(channel_slice_outside)
% 
%                         axis image
%                         set(gca,'Ydir','Reverse')
% 
%                         set(gca,'Color','white')
%                         set(gca,'XTick',[])
%                         set(gca,'YTick',[])
%                 
%                 suptitle(['p value (sig over background) = ' num2str(p_value)])
%                 colormap gray
%                 
%                 pause
%             end
%             
%             figure(13423)
%             clf
%             
%                 subtightplot(1,3,1)
%                     hold all
%                     
%                     imagesc(cur_slice(:,:,options.nuclei_channel))
%                     
%                     axis image
%                     set(gca,'Ydir','Reverse')
% 
%                     set(gca,'Color','white')
%                     set(gca,'XTick',[])
%                     set(gca,'YTick',[])
%                     
%                 subtightplot(1,3,2)
%                     hold all
%                     
%                     imagesc(cur_slice(:,:,options.cell_channels(1)))
%                     
%                     axis image
%                     set(gca,'Ydir','Reverse')
% 
%                     set(gca,'Color','white')
%                     set(gca,'XTick',[])
%                     set(gca,'YTick',[])
%                     
%                 subtightplot(1,3,3)
%                     hold all
%                     
%                     imagesc(cur_slice(:,:,options.cell_channels(2)))
%                     
%                     axis image
%                     set(gca,'Ydir','Reverse')
% 
%                     set(gca,'Color','white')
%                     set(gca,'XTick',[])
%                     set(gca,'YTick',[])
%                     
%             colormap gray
%             
%             set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
%             set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
%             set(gcf, 'color', 'white');
%             
%             cell_segmentation_results_struct_out = cell_segmentation_results_struct;
%             
%             pause
        end
    end
    
end

