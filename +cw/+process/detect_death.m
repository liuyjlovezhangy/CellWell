function cell_death_results = detect_death( well_tracking_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, options )
    propts = options.processing_options;
    
    num_wells = numel(well_tracking_results_struct.wells);
    
    for well_idx = 1:num_wells
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        
        for channel_idx = options.cell_channels
            
            channel_cells = cell_tracking_results_struct.linked_object_cells{well_idx,channel_idx};
            
            for cell_idx = 1:numel(channel_cells)
                cell = channel_cells(cell_idx);
                
                for frame_idx = 1:cell.end_frame - cell.start_frame + 1
                
                    props = cell.cell_props(frame_idx);

                    start_i = ceil(props.BoundingBox(2));
                    end_i = ceil(start_i + props.BoundingBox(4)-1);

                    start_j = ceil(props.BoundingBox(1));
                    end_j = ceil(start_j + props.BoundingBox(3)-1);

                    cell_image = props.FilledImage;

                    cur_slice = cur_well_im(:,:,cell.start_frame + frame_idx - 1,options.signal_channels);

%                     cur_slice = wiener2(cur_slice);
                    
                    sig_slice = cur_slice(start_i:end_i,start_j:end_j);
                    sig_slice_thresh = sig_slice;
                    sig_slice_thresh(cell_image == 0) = 0;
                    
                    cur_nuclei_mask = cell_segmentation_results_struct.cell_masks{well_idx}(:,:,cell.start_frame + frame_idx - 1,options.nuclei_channel);
                    
                    for_slice = cur_slice;
                    for_slice(cur_nuclei_mask(:) == 0) = 0;
                    for_slice = reshape(for_slice,size(cur_nuclei_mask));
                    
                    back_slice = cur_slice;
                    back_slice(cur_nuclei_mask(:) ~= 0) = 0;
                    back_slice = reshape(back_slice,size(cur_nuclei_mask));

                    % calculate Mann-Whitney of signal compared to background
                    
                    p_value = ranksum(sig_slice_thresh(sig_slice_thresh ~= 0),back_slice(back_slice ~= 0));
                    
                    if p_value < 1e-12
                    
                        figure(82493)
                        clf
                        hist(sig_slice_thresh(sig_slice_thresh ~= 0),100)
                        xlim([0 0.2])
                        
                        figure(82393)
                        clf
                        hist(back_slice(back_slice ~= 0),100)
                        xlim([0 0.2])
                        
                        figure(13742)
                        clf
                            subtightplot(2,4,1)
                                hold all

                                imagesc(cur_slice)

                                axis image
                                set(gca,'Ydir','Reverse')

                                set(gca,'Color','white')
                                set(gca,'XTick',[])
                                set(gca,'YTick',[])
                                
                            subtightplot(2,4,2)
                                hold all

                                imagesc(for_slice)

                                axis image
                                set(gca,'Ydir','Reverse')

                                set(gca,'Color','white')
                                set(gca,'XTick',[])
                                set(gca,'YTick',[])

                            subtightplot(2,4,3)
                                hold all

                                imagesc(sig_slice_thresh)

                                axis image
                                set(gca,'Ydir','Reverse')

                                set(gca,'Color','white')
                                set(gca,'XTick',[])
                                set(gca,'YTick',[])
                                
                            subtightplot(2,4,4)
                                hold all

                                imagesc(back_slice)

                                axis image
                                set(gca,'Ydir','Reverse')

                                set(gca,'Color','white')
                                set(gca,'XTick',[])
                                set(gca,'YTick',[])
                                
                            subplot(2,2,3)
                            
                                hist(sig_slice_thresh(:))
                                
                                xlim([0 1])
                                
                                title('Signal histogram')
                                ylabel('Count')
                                xlabel('Intensity')
                                
                            subplot(2,2,4)
                            
                                hist(back_slice(:))
                                
                                xlim([0 1])
                                
                                title('Background histogram')
                                ylabel('Count')
                                xlabel('Intensity')
                                

                        colormap gray

                        set(gcf,'color','w')
                        suptitle(['P-value: ' num2str(p_value)])
%                         suptitle([num2str(well_idx) ' ' num2str(channel_idx) ' ' num2str(cell_idx) ' ' num2str(p_value)])
                        
                        pause
                    end
                end
            end
        end
    end
end