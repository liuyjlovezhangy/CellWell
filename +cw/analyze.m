function analyze(options)
    disp('Analysis running...')

    %%% Load movie
    
    full_filename = options.filename;
    
    disp('Loading movie...')
    
    im = importdata([full_filename '__analysis_results/input_movie.mat']);
    im = mat4D_to_gray(im);
    
    well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    well_tracking_results_struct = cw.process.track_wells( im, well_segmentation_results_struct, options );
    signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
    cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
    cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    cell_interaction_results_struct = importdata([full_filename '__analysis_results/cell_interactions.mat']);
    
    interaction_results(cell_interaction_results_struct,options);
    
%     cell_location_distributions(cell_tracking_results_struct,options);
end

function interaction_results(cell_interaction_results_struct,options)
    num_wells = numel(cell_interaction_results_struct.interactions);
    
    warning('interaction_results: hack where it assumes 2 cell channels')
    
    unique_interaction_counts = zeros(num_wells, 3);
    
    % wells where there are no <x> cells should not contribute to histogram
    
    
    
    for well_idx = 1:num_wells
        well_interactions = cell_interaction_results_struct.interactions{well_idx};
        
        for interaction_idx = 1:numel(well_interactions)
            int = well_interactions(interaction_idx);
            
            % how many unique interactions for these two cells
            
            changes = int.interacting_flag(2:end) - int.interacting_flag(1:end-1);
            UIC = sum(changes == -1) + 1;
            
            channel_offset = options.cell_channels(1) - 1;
            
            if int.channel1 == int.channel2
                UIC_col = int.channel1 - channel_offset;
            else
                UIC_col = 3;
            end
            
            unique_interaction_counts(well_idx, UIC_col) = unique_interaction_counts(well_idx, UIC_col) + UIC;
            
            % how long did each interaction last?
            
            
        end

    end
    
    UIC_1_hist_loc = 0:max(unique_interaction_counts(:,1))+1;
    UIC_2_hist_loc = 0:max(unique_interaction_counts(:,2))+1;
    UIC_combo_hist_loc = 0:max(unique_interaction_counts(:,3))+1;

    UIC_1_hist = histc(unique_interaction_counts(:,1),UIC_1_hist_loc);
    UIC_2_hist = histc(unique_interaction_counts(:,2),UIC_2_hist_loc);
    UIC_combo_hist = histc(unique_interaction_counts(:,3),UIC_combo_hist_loc);
    
    figure(134287)
    clf
        subplot(1,3,1)
            hold all

            plot(UIC_1_hist_loc,UIC_1_hist,'.-','Color','r','LineWidth',3,'MarkerSize',40)
            plot(UIC_2_hist_loc,UIC_2_hist,'.-','Color','g','LineWidth',3,'MarkerSize',40)
            plot(UIC_combo_hist_loc,UIC_combo_hist,'.-','Color','b','LineWidth',3,'MarkerSize',40)

            xlabel('Number of unique interactions')
            ylabel('# wells')
            title('Number of unique cell-cell interactions found per well')

            box on
            grid on
            
        subplot(1,3,3)
            hold all
            
            plot(-1,-1,'.-','Color','r','LineWidth',3,'MarkerSize',40)
            plot(-1,-1,'.-','Color','g','LineWidth',3,'MarkerSize',40)
            plot(-1,-1,'.-','Color','b','LineWidth',3,'MarkerSize',40)
            
            xlim([0 1])
            ylim([0 1])
            axis off
            
            comb_string = [options.channel_labels{options.cell_channels(1)} ' <-> ' options.channel_labels{options.cell_channels(2)} ' interaction'];
            
            legend([options.channel_labels(options.cell_channels), {comb_string}])
        
%     suptitle('Cell-cell interaction statistics')
            
    set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
    set(gcf, 'color', 'white')
end

function cell_location_distributions(cell_tracking_results_struct,options)
    warning('cell_location_distributions should really be made colorblind-compatible...')
    warning('cell_location_distributions assumes cell channels are continuous at the moment...')
    
    cell_channels = sort(options.cell_channels);
    
    ctrs = cell_tracking_results_struct;
    
    num_wells = size(ctrs.cell_tracks,1);

    well_center = [options.well_width/2; options.well_height/2];
    max_dist = eucl_dist([0;0], well_center);
    
    cell_channel_colors = {'r','g'};
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Simulate the positions of uniform random samples in a well
    
    random_localization_count = 1e5;
    
    %%% Everywhere
    
    random_x_loc = options.well_width * rand(1,random_localization_count);
    random_y_loc = options.well_height * rand(1,random_localization_count);
    
    random_distances = eucl_dist([random_x_loc;random_y_loc], repmat(well_center,[1,random_localization_count]));
    random_angles = adjust_atan2(random_y_loc - well_center(2),random_x_loc - well_center(1)) * 360 / (2*pi);
    
    [random_distances_hist,random_distances_hist_loc] = histogram(random_distances);
    random_distances_hist = mat2gray(random_distances_hist);
    
    [random_angles_hist,random_angles_hist_loc] = histogram(random_angles);
    random_angles_hist = mat2gray(random_angles_hist);
    
    %%% Randomly attached to an edge
    
    rxle1 = options.well_width * rand(1,random_localization_count);
    rxle2 = options.well_width * rand(1,random_localization_count);
    rxle3 = ones(1,random_localization_count);
    rxle4 = options.well_width * ones(1,random_localization_count) - 1;
    
    ryle1 = options.well_height * ones(1,random_localization_count) - 1;
    ryle2 = ones(1,random_localization_count);
    ryle3 = options.well_height * rand(1,random_localization_count);
    ryle4 = options.well_height * rand(1,random_localization_count);
    
    random_x_loc_edge = [rxle1, rxle2, rxle3, rxle4];              
    random_y_loc_edge = [ryle1, ryle2, ryle3, ryle4];
    
    random_distances_edge = eucl_dist([random_x_loc_edge;random_y_loc_edge], repmat(well_center,[1,numel(random_x_loc_edge)]));
    random_angles_edge = adjust_atan2(random_y_loc_edge - well_center(2),random_x_loc_edge - well_center(1)) * 360 / (2*pi);
    
    [random_distances_edge_hist,random_distances_edge_hist_loc] = histogram(random_distances_edge);
    random_distances_edge_hist = mat2gray(random_distances_edge_hist);
    
    [random_angles_edge_hist,random_angles_edge_hist_loc] = histogram(random_angles_edge);
	random_angles_edge_hist = mat2gray(random_angles_edge_hist);
    
    %%% Compile cell localizations across each well & channel
    
    channel_localizations = cell(1,numel(cell_channels));
    channel_angles = channel_localizations;
    
    channel_distances = channel_localizations;
    channel_distances_hist = channel_localizations;
    channel_distances_hist_loc = channel_localizations;
    
    channel_angles_hist = channel_localizations;
    channel_angles_hist_loc = channel_localizations;
    
    for channel_idx = cell_channels
        for well_idx = 1:num_wells
            temp_tracks = ctrs.cell_tracks{well_idx,channel_idx};
            temp_localizations = [];
            
            for track_idx = 1:numel(temp_tracks)
                track = temp_tracks{track_idx};
                
                cur_localizations = track(:,~isnan(sum(track,1)));
                
                temp_localizations = [temp_localizations, cur_localizations];
            end
            
            if ~isempty(temp_localizations)
                channel_localizations{channel_idx - cell_channels(1) + 1} = [channel_localizations{channel_idx - cell_channels(1) + 1}, temp_localizations];
                
                temp_distances = eucl_dist(temp_localizations, repmat(well_center,[1,size(temp_localizations,2)]));
            
                temp_angles = adjust_atan2(temp_localizations(2,:) - well_center(2),temp_localizations(1,:) - well_center(1))  * 360 / (2*pi);
                
                channel_distances{channel_idx - cell_channels(1) + 1} = [channel_distances{channel_idx - cell_channels(1) + 1}, temp_distances];
                channel_angles{channel_idx - cell_channels(1) + 1} = [channel_angles{channel_idx - cell_channels(1) + 1}, temp_angles];
            end
        end

        [channel_distances_hist{channel_idx - cell_channels(1) + 1},...
            channel_distances_hist_loc{channel_idx - cell_channels(1) + 1}] = ...
            histogram(channel_distances{channel_idx - cell_channels(1) + 1});
        channel_distances_hist{channel_idx - cell_channels(1) + 1} = ...
            mat2gray(channel_distances_hist{channel_idx - cell_channels(1) + 1});
        
        [channel_angles_hist{channel_idx - cell_channels(1) + 1},...
            channel_angles_hist_loc{channel_idx - cell_channels(1) + 1}] = ...
            histogram(channel_angles{channel_idx - cell_channels(1) + 1});
        channel_angles_hist{channel_idx - cell_channels(1) + 1} = ...
            mat2gray(channel_angles_hist{channel_idx - cell_channels(1) + 1});
    end
    
    %%% Plot
    
    figure(13542)
    clf
    
        subtightplot(2,3,1,0.05)
            hold all

            plot(random_x_loc(1:1e4),random_y_loc(1:1e4),'.','Color',[0.5 0.5 0.5],'MarkerSize',4)
            
            for channel_idx = 1:numel(cell_channels)
                plot(channel_localizations{channel_idx}(1,:),channel_localizations{channel_idx}(2,:),'.','Color',cell_channel_colors{channel_idx},'MarkerSize',10)
            end
            
            plot([well_center(1),well_center(1)],[0,options.well_height],'--w','LineWidth',3)
            plot([0,options.well_width],[well_center(2),well_center(2)],'--w','LineWidth',3)
            plot(well_center(1),well_center(2),'ow','MarkerSize',20,'LineWidth',5)
            
            circle_segment_arrow(well_center(1)+5,well_center(2),10,90,50,'w');
            
            xlim([0 options.well_width])
            ylim([0 options.well_height])

            title('Localizations of cells across all wells')
            
            axis image
            set(gca,'Ydir','Reverse')
            set(gca, 'color', 'black')
            
            box on
            
        subtightplot(2,3,2,0.05)
            hold all

            for channel_idx = 1:numel(cell_channels)
                plot(channel_distances_hist_loc{channel_idx},channel_distances_hist{channel_idx},'-','Color',cell_channel_colors{channel_idx},'LineWidth',3)
            end
            
            plot(random_distances_hist_loc,random_distances_hist,'--','Color',[0.5 0.5 0.5],'LineWidth',3)

            xlim([0 max_dist])
            ylim([0 1])%max([channel_distances_hist{1}, channel_distances_hist{2}, random_distances_hist, random_distances_edge_hist])])
            
            ylabel('Normalized density')
            xlabel('Distance from center')
            title('Distribution of distances of cells from center')

            box on
%             grid on
            
        subtightplot(2,3,3,0.05)
            hold all
            
            for channel_idx = 1:numel(cell_channels)
                plot(channel_angles_hist_loc{channel_idx},channel_angles_hist{channel_idx},'-','Color',cell_channel_colors{channel_idx},'LineWidth',3)
            end
            
            plot(random_angles_hist_loc,random_angles_hist,'--','Color',[0.5 0.5 0.5],'LineWidth',3)
%             plot(random_angles_edge_hist_loc,random_angles_edge_hist,'-b','LineWidth',3)
            
            xlim([0 360])
            ylim([0 1])
            
            ylabel('Normalized density')
            xlabel('Angle from horizontal (degrees)')
            title('Distribution of angle from horizontal')
            
            set(gca,'xtick',0:45:360)
            
            box on
%             grid on
            
        subtightplot(2,3,4,0.05)
            hold all
            
            for channel_idx = 1:numel(cell_channels)
                plot(-1,-1,'-','LineWidth',3,'Color',cell_channel_colors{channel_idx})
            end
            
            plot(-1,-1,'--','Color',[0.5 0.5 0.5],'LineWidth',3)
            
            legend([options.channel_labels(cell_channels), 'Uniform distribution'])
            
            axis off
        
    set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
    set(gcf, 'color', 'white')
end

function circle_segment_arrow(x,y,r,degrees,nsegments,color)  
    if nargin<4
        nsegments=50;
    end
    
    th = linspace(0,2*pi*degrees/360,nsegments);
    xunit = r * cos(th) + x;
    yunit = r * sin(th) + y;
    plot(xunit(1:end-20), yunit(1:end-20),'LineWidth',6,'Color',color);
    
    plot_arrow(xunit(end-20),yunit(end-20),xunit(end),yunit(end),'linewidth',4,'color',color,'facecolor',color,'headwidth',1000,'headheight',1000,'edgecolor',color);
    
%     set(gca,'children',flipud(get(gca,'children')))
end