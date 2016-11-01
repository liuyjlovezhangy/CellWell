function analyze(options)
    disp('Analysis running...')

    %%% Load movie
    
    full_filename = options.filename;
    
    disp('Loading movie...')
    
    im = importdata([full_filename '__analysis_results/input_movie.mat']);
    im = mat4D_to_gray(im);
    
    well_segmentation_results_struct = [];
    well_tracking_results_struct = [];
    signal_detection_results_struct = [];
    cell_segmentation_results_struct = [];
    cell_tracking_results_struct = [];
    
    disp('Loading data...')
    
%     if isempty(well_segmentation_results_struct)
%         well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
%     end
%         
%     if isempty(well_tracking_results_struct)
%        well_tracking_results_struct = importdata([full_filename '__analysis_results/well_tracking.mat']);
%     end
%     
%     if isempty(signal_detection_results_struct)
%        signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
%     end
%     
%     if isempty(cell_segmentation_results_struct)
%        cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
%     end
    
    if isempty(cell_tracking_results_struct)
       cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    end
    
    cell_location_distributions(cell_tracking_results_struct,options);
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