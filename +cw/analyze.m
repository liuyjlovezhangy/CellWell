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
    well_center = [options.well_width; options.well_height];
    max_dist = eucl_dist([0;0], well_center);
    
    %%% Simulate the positions of uniform random samples in a well
    
    random_localization_count = 1e5;
    
    random_x_loc = options.well_width * rand(1,random_localization_count);
    random_y_loc = options.well_height * rand(1,random_localization_count);
    
    random_distances = eucl_dist([random_x_loc;random_y_loc], repmat(well_center,[1,random_localization_count]));
    [random_distances_hist,random_distances_hist_loc] = histogram(random_distances);
    random_distances_hist = random_distances_hist / sum(random_distances_hist);
    
    %%% Compile cell localizations across each well & channel
    
    channel_localizations = cell(1,numel(options.cell_channels));
    
    %%% Plot
    
    figure(13542)
    clf
    
        subplot(1,2,1)
            hold all

            plot(random_x_loc(1:1e4),random_y_loc(1:1e4),'.b','MarkerSize',10,'LineWidth',3)

            xlim([0 options.well_width])
            ylim([0 options.well_height])

            box on
            
        subplot(1,2,2)
            hold all

            plot(random_distances_hist_loc,random_distances_hist,'-b','LineWidth',3)

%             xlim([0 options.well_width])
%             ylim([0 options.well_height])

%             line([max_dist, max_dist], ylim, 'Color','k','LineWidth',3)


        
    set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
    set(gcf, 'color', 'white');
end