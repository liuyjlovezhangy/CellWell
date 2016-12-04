function cell_interaction_results_struct = detect_interactions( well_tracking_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, options )
    
    num_wells = size(cell_tracking_results_struct.cell_tracks,1);
    num_frames = size(well_tracking_results_struct.well_tracks_mat_shifted,2);
    
    cell_interaction_results_struct = [];
    cell_interaction_results_struct.interactions = cell(1,num_wells);
    
    %%% sweep through tracks and find other tracks within a certain distance

    for well_idx = 1:num_wells

        cell_mask = cell_segmentation_results_struct.cell_masks{well_idx};
        
        well_interactions = [];
        
        for source_channel_idx = options.cell_channels
            source_tracks_cell = cell_tracking_results_struct.cell_tracks{well_idx,source_channel_idx};
            
            for target_channel_idx = options.cell_channels(options.cell_channels >= source_channel_idx)
                target_tracks_cell = cell_tracking_results_struct.cell_tracks{well_idx,target_channel_idx};
                
                % go through each track in the source cell and compare to
                % all tracks in the target cell
                
                for source_track_idx = 1:numel(source_tracks_cell)
                    source_track = source_tracks_cell{source_track_idx};
                    
                    if source_channel_idx == target_channel_idx
                        target_vec = source_track_idx+1:numel(target_tracks_cell);
                    else
                        target_vec = 1:numel(target_tracks_cell);
                    end
                        
                    for target_track_idx = target_vec
                        target_track = target_tracks_cell{target_track_idx};
                        
                        cur_objects_sc = cell_tracking_results_struct.linked_object_cells{well_idx,source_channel_idx};
                        cur_objects_tc = cell_tracking_results_struct.linked_object_cells{well_idx,target_channel_idx};
                        
                        cur_source_object_props = cur_objects_sc(source_track_idx).cell_props;
                        cur_target_object_props = cur_objects_tc(target_track_idx).cell_props;
                        
                        source_radii = NaN * ones(1,num_frames);%[cur_source_object_props.EquivDiameter]/2
                        target_radii = NaN * ones(1,num_frames);%[cur_target_object_props.EquivDiameter]/2
                        
                        source_radii(cur_objects_sc(source_track_idx).start_frame:cur_objects_sc(source_track_idx).end_frame) = [cur_source_object_props.EquivDiameter]/2;
                        target_radii(cur_objects_tc(target_track_idx).start_frame:cur_objects_tc(target_track_idx).end_frame) = [cur_target_object_props.EquivDiameter]/2;
                        
                        max_dists = sum([source_radii;target_radii],1) + options.interaction_tolerance;
                        
                        distances = eucl_dist(source_track,target_track);
                        interacting_flag = distances <= max_dists;

                        if all(isnan(distances))
                            continue
                        end

                            interaction.distances = distances;
                            interaction.interacting_flag = interacting_flag;
                            interaction.id1 = source_track_idx;
                            interaction.id2 = target_track_idx;
                            interaction.track1 = source_track;
                            interaction.track2 = target_track;
                            interaction.channel1 = source_channel_idx;
                            interaction.channel2 = target_channel_idx;
                        
                            well_interactions = [well_interactions, interaction];

                        if 0
                        
                            for frame_idx = 1:num_frames

                                figure(10394)
                                clf

                                    subtightplot(2,3,1)
                                        hold all

                                        imagesc(cell_mask(:,:,frame_idx,source_channel_idx))

                                        plot(source_track(1,frame_idx),source_track(2,frame_idx),'mo','MarkerSize',15,'LineWidth',3)
                                        plot(target_track(1,frame_idx),target_track(2,frame_idx),'mx','MarkerSize',15,'LineWidth',3)

                                        if interacting_flag(frame_idx)
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'g--','LineWidth',4)
                                        else
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'r--','LineWidth',4)
                                        end

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                        title(['Source channel: ' num2str(source_channel_idx) ])

                                    subtightplot(2,3,4)
                                        hold all

                                        imagesc(well_tracking_results_struct.wells(well_idx).im_well(:,:,frame_idx,source_channel_idx))

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                    subtightplot(2,3,2)
                                        hold all

                                        imagesc(cell_mask(:,:,frame_idx,target_channel_idx))

                                        plot(source_track(1,frame_idx),source_track(2,frame_idx),'mo','MarkerSize',15,'LineWidth',3)
                                        plot(target_track(1,frame_idx),target_track(2,frame_idx),'mx','MarkerSize',15,'LineWidth',3)

                                        if interacting_flag(frame_idx)
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'g--','LineWidth',4)
                                        else
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'r--','LineWidth',4)
                                        end

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                        title(['Target channel: ' num2str(target_channel_idx) ])

                                    subtightplot(2,3,5)
                                        hold all

                                        imagesc(well_tracking_results_struct.wells(well_idx).im_well(:,:,frame_idx,target_channel_idx))

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                    subtightplot(2,3,3)
                                        hold all

                                        imagesc(cell_mask(:,:,frame_idx,options.nuclei_channel))

                                        plot(source_track(1,frame_idx),source_track(2,frame_idx),'mo','MarkerSize',15,'LineWidth',3)
                                        plot(target_track(1,frame_idx),target_track(2,frame_idx),'mx','MarkerSize',15,'LineWidth',3)

                                        if interacting_flag(frame_idx)
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'g--','LineWidth',4)
                                        else
                                            plot([source_track(1,frame_idx),target_track(1,frame_idx)],[source_track(2,frame_idx),target_track(2,frame_idx)],...
                                                'r--','LineWidth',4)
                                        end

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                        title('Nuclei channel')

                                    subtightplot(2,3,6)
                                        hold all

                                        imagesc(well_tracking_results_struct.wells(well_idx).im_well(:,:,frame_idx,options.nuclei_channel))

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

                                        freezeColors

                                    suptitle(['Well: ' num2str(well_idx) ' source channel: ' num2str(source_channel_idx) ' source cell: ' num2str(source_track_idx) ' frame: ' ...
                                        num2str(frame_idx) ' target channel: ' num2str(target_channel_idx), ' target cell: ' num2str(target_track_idx)])

                                    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
                                    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
                                    set(gcf, 'color', 'white');
                                    drawnow
                                pause
                            end
                        end
                    end
                end
            end
        end
        
        
        cell_interaction_results_struct.interactions{well_idx} = well_interactions;
    end

end