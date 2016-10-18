function final_analysis(well_tracking_results_struct,signal_detection_results_struct,cell_segmentation_results_struct,cell_tracking_results_struct)

    inm = signal_detection_results_struct.is_noise_matrix;

    num_wells = size(inm,1);
    num_frames = size(inm,3);
    
    num_total_channels = size(inm,2) + 1;
    num_cell_channels = 2;
    num_signal_channels = num_total_channels - num_cell_channels;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% How many wells had no signal?
    
    no_signal_percentages_wells = (num_wells - sum(any(~inm,3),1)) / num_wells * 100;
    
    no_signal_percentages_frames_dist = sum(inm,3) / num_frames * 100;
    
    no_signal_percentages_frames_means = mean(no_signal_percentages_frames_dist,1);
    no_signal_percentages_frames_stds = std(no_signal_percentages_frames_dist,0,1);
    
    datasets = {};
    
    for set_idx = 1:num_total_channels-1
        datasets{set_idx} = no_signal_percentages_frames_dist(:,set_idx);
    end
    
    [NSPF_p_values,NSPF_groups] = permutive_ttesting( datasets );
    
    figure(134234)
    clf
        % bar and error bar alignment
        
        num_bars = numel(no_signal_percentages_wells);
    
        subplot(1,2,1)
            hold all

            bar(1, no_signal_percentages_wells(1), 'facecolor', 'r'); 
            bar(2, no_signal_percentages_wells(2), 'facecolor', 'g'); 
            bar(3, no_signal_percentages_wells(3), 'facecolor', 'b');

            title('Percentage of wells with no signal')
            ylabel('Percentage (%)')
            xlabel('Imaging channel')

            ylim([0 100])

            set(gca,'XTick',1:num_bars)
            set(gca,'XTickLabel',{'K562','CD19_CART','SYTOX'})
            
        subplot(1,2,2)
            hold all

            bar(1, no_signal_percentages_frames_means(1), 'facecolor', 'r'); 
            bar(2, no_signal_percentages_frames_means(2), 'facecolor', 'g'); 
            bar(3, no_signal_percentages_frames_means(3), 'facecolor', 'b');
            
            errorbar(1:num_bars, no_signal_percentages_frames_means, no_signal_percentages_frames_stds, 'k', 'LineWidth',5, 'linestyle', 'none');
            
            sigstar(NSPF_groups,NSPF_p_values);
            
            title('Percentage of frames with no signal')
            ylabel('Percentage (%)')
            xlabel('Imaging channel')
            
            yl = ylim;
            ylim([0,yl(2)])
            set(gca,'YTick',0:10:100)

%             ylim([0 100])

            set(gca,'XTick',1:num_bars)
            set(gca,'XTickLabel',{'K562','CD19_CART','SYTOX'})
        
%     suptitle('Quantification of numbers of wells and frames with no signal')
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');

    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% How many cells of each type max per well?
    
    max_number_objects_dist = zeros(num_wells,num_cell_channels);
    
    for well_idx = 1:num_wells
        for frame_idx = 1:num_frames
            for channel_idx = 1:num_cell_channels
                current_objects = cell_segmentation_results_struct.detected_cell_props_final{well_idx,frame_idx,channel_idx};
                
                max_number_objects_dist(well_idx,channel_idx) = max([max_number_objects_dist(well_idx,channel_idx), numel(current_objects)]);
            end
        end
    end
    
    max_number_objects_mean = mean(max_number_objects_dist,1);
    max_number_objects_std = std(max_number_objects_dist,0,1);
    
    datasets = {};
    
    for set_idx = 1:num_cell_channels
        datasets{set_idx} = max_number_objects_dist(:,set_idx);
    end
    
    [MNO_p_values,MNO_groups] = permutive_ttesting( datasets );
    
    num_bars = num_cell_channels;
    
    figure(135423)
    clf
    
        subplot(1,1,1)
            hold all

            bar(1, max_number_objects_mean(1), 'facecolor', 'r'); 
            bar(2, max_number_objects_mean(2), 'facecolor', 'g'); 
            
            errorbar(1:num_bars, max_number_objects_mean, max_number_objects_std, 'k', 'LineWidth',5, 'linestyle', 'none');
            
            sigstar(MNO_groups,MNO_p_values);
            
            title('Number of cells in a well')
            ylabel('Number of cells')
            xlabel('Imaging channel')
            
            yl = ylim;
            ylim([0,yl(2)])

            set(gca,'XTick',1:num_bars)
            set(gca,'XTickLabel',{'K562','CD19_CART'})
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');

    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Trajectory statistics
    
    cell_tracks = cell_tracking_results_struct.cell_tracks;
    
    track_lengths = cell(1,num_cell_channels);
    
    for well_idx = 1:num_wells
        for channel_idx = 1:num_cell_channels
            cur_tracks = cell_tracks{well_idx,channel_idx};
            
            cur_tlens = [];
            
            for track_idx = 1:numel(cur_tracks)
                track = cur_tracks{track_idx};
                num_points = numel(find(~isnan(sum(track,1))));
                
                cur_tlens = [cur_tlens, num_points];
            end
            
            track_lengths{channel_idx} = [track_lengths{channel_idx}, cur_tlens];
        end
    end
    
    track_len_means = zeros(1,num_cell_channels);
    track_len_stds = zeros(1,num_cell_channels);
    
    for channel_idx = 1:num_cell_channels
        track_len_means(channel_idx) = mean(track_lengths{channel_idx});
        track_len_stds(channel_idx) = std(track_lengths{channel_idx});
    end
    
    [TL_p_values,TL_groups] = permutive_ttesting( track_lengths );
    
    num_bars = num_cell_channels;
    
    figure(113923)
    clf
    
        subplot(1,1,1)
            hold all

            bar(1, track_len_means(1), 'facecolor', 'r'); 
            bar(2, track_len_means(2), 'facecolor', 'g'); 
            
            errorbar(1:num_bars, track_len_means, track_len_stds, 'k', 'LineWidth',5, 'linestyle', 'none');
            
            sigstar(TL_groups,TL_p_values);
            
            title('Length of trajectories')
            ylabel('Track length [frames]')
            xlabel('Imaging channel')
            
            yl = ylim;
            ylim([0,yl(2)])

            set(gca,'XTick',1:num_bars)
            set(gca,'XTickLabel',{'K562','CD19_CART'})
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');

    drawnow
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Cell object statistics
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Stats about wall-hugging cells
    
    
end