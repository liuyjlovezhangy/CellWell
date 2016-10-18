function test_segmentation_clean_up

    % load test data
    
    well_tracking_results_struct = importdata('code/test_cases/well_tracking_test_data.mat');
    cell_segmentation_results_struct = importdata('code/test_cases/cell_segmentation_test_data.mat');
    
    % run through the watershedded and non-watershedded object maps to pick
    % out objects that were removed inappropriately
    
    cell_segmentation_results_corrected_struct = clean_up_cell_segmentation( cell_segmentation_results_struct, [] );
    
    save('code/test_cases/cell_segmentation_corrected_test_data.mat','cell_segmentation_results_corrected_struct');
    
    % plot original and corrected results
    
    wells = well_tracking_results_struct.wells;
    cur_well_img = wells(1).im_well;
    
    nowater_masks = cell_segmentation_results_corrected_struct.cell_masks_nowater;
    water_masks = cell_segmentation_results_corrected_struct.cell_masks_water;
    final_masks = cell_segmentation_results_corrected_struct.cell_masks_final;
    
    nowater_cells = cell_segmentation_results_corrected_struct.detected_cell_props_nowater;
    water_cells = cell_segmentation_results_corrected_struct.detected_cell_props_water;
    final_cells = cell_segmentation_results_corrected_struct.detected_cell_props_final;
    
    num_frames = size(nowater_cells,2);
    num_channels = size(nowater_cells,3);
    
    num_colors = 2^16;
    
    color_space = linspace(0, 1, num_colors)';
    z = zeros(num_colors,1);
    
    r_map = [color_space, z, z];
    g_map = [z, color_space, z];
    b_map = [z, z, color_space];
    
    cmaps = {r_map,g_map,b_map};
    
    for frame_idx = 1:num_frames
        figure(13423)
        clf
        
        for channel_idx = 1:num_channels
            water_mask_px = water_masks{1}(:,:,frame_idx,channel_idx);
            
            if all(water_mask_px(:) == 0)
                zeromap_water = [0 0 1];
            else
                zeromap_water = [.5 .5 .5];
            end
            
            final_mask_px = final_masks{1}(:,:,frame_idx,channel_idx);
            
            if all(final_mask_px(:) == 0)
                zeromap_final = [0 0 1];
            else
                zeromap_final = [.5 .5 .5];
            end
            
            subtightplot(num_channels,4,4*(channel_idx-1)+1)
                hold all
                
                imagesc(cur_well_img(:,:,frame_idx,channel_idx))
                
                colormap(cmaps{channel_idx})
                
                freezeColors
                
                axis image
                axis tight
                set(gca,'Ydir','Reverse')
                axis off
                
                if channel_idx == 1
                    title('Original image')
                end
                ylabel(['Channel: ' num2str(channel_idx)])
            
            subtightplot(num_channels,4,4*(channel_idx-1)+2)
                hold all
                
                nowater_slice = nowater_masks{1}(:,:,frame_idx,channel_idx);
                if all(nowater_slice(:)==0)
                    slice = zeros(size(nowater_slice,1),size(nowater_slice,2),3);
                    slice(:,:,3) = 1;
                    imshow(slice)
                else
                    rgb = label2rgb(nowater_slice+1,'jet');
                    imshow(rgb)
                end
                
                axis image
                axis tight
                set(gca,'Ydir','Reverse')
                axis off
                
                if channel_idx == 1
                    title('Non-watershed segmentation')
                end
                
            subtightplot(num_channels,4,4*(channel_idx-1)+3)
                hold all
                
                rgb = label2rgb(water_masks{1}(:,:,frame_idx,channel_idx),'jet',zeromap_water);
                imshow(rgb)
                
                axis image
                axis tight
                set(gca,'Ydir','Reverse')
                axis off
                
                if channel_idx == 1
                    title('Watershed segmentation')
                end
                
            subtightplot(num_channels,4,4*(channel_idx-1)+4)
                hold all
                
                rgb = label2rgb(final_masks{1}(:,:,frame_idx,channel_idx),'jet',zeromap_final);
                imshow(rgb)
                
                axis image
                axis tight
                set(gca,'Ydir','Reverse')
                axis off
                
                if channel_idx == 1
                    title('Final segmentation')
                end
        end
        
        suptitle(['Cell segmentation correction test. Frame: ' num2str(frame_idx) ' of ' num2str(num_frames)])

        set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');

        drawnow
    end

end
