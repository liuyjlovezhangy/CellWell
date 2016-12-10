function analyze(options)
    disp('Analysis running...')

    %%% Load movie
    
    full_filename = options.filename;
    
    disp('Loading movie...')
    
    im = importdata([full_filename '__analysis_results/input_movie.mat']);
    im = mat4D_to_gray(im);
    
    disp matfiles
    
    well_segmentation_results_struct = importdata([full_filename '__analysis_results/well_segmentation.mat']);
    well_tracking_results_struct = cw.process.track_wells( im, well_segmentation_results_struct, options );
    signal_detection_results_struct = importdata([full_filename '__analysis_results/noise_detection.mat']);
    cell_segmentation_results_struct = importdata([full_filename '__analysis_results/cell_segmentation.mat']);
    cell_tracking_results_struct = importdata([full_filename '__analysis_results/cell_tracking.mat']);
    cell_interaction_results_struct = importdata([full_filename '__analysis_results/cell_interactions.mat']);
    
    disp done
    
    track_cdfs( cell_tracking_results_struct,options );
    
%     interaction_results(cell_tracking_results_struct, cell_interaction_results_struct,options);
    
%     cell_location_distributions(cell_tracking_results_struct,options);
end

function track_cdfs( cell_tracking_results_struct,options )

    F_test_alpha = 0.01;

    delta_t = options.time_step;
    umperpx = options.pixel_size;

    num_iterations = 30;
    max_states = 2;
    
    ct = cell_tracking_results_struct.cell_tracks;
    
    for channel_idx = options.cell_channels
        
        tracks = {};

        for well_idx = 1:size(ct,1)
            tracks = [tracks, ct{well_idx,channel_idx}];
        end

        %%%%%%%%%%%%%%%%%%

        D_cell = cell(1,max_states);
        weight_cell = cell(1,max_states);

        steps_cell = {};
        all_delta_rs = [];

        for track_idx = 1:numel(tracks)
            track = tracks{track_idx};

            steps = track(:,2:end) - track(:,1:end-1);
            
%             steps = 1*randn(2,20);

            steps_cell = [steps_cell, steps];

            delta_rs = sqrt( steps(1,:) .^ 2 + steps(2,:) .^2 );

            delta_rs = delta_rs(~isnan(delta_rs));

            all_delta_rs = [all_delta_rs, delta_rs];
        end

        all_delta_rs = all_delta_rs * umperpx;

%         [data_cdf,x_locations] = zecdf(all_delta_rs);
        
        [data_cdf,x_locations] = ecdf(all_delta_rs);
        
%         [data_cdf data_cdf2]
        
%         return

        %%%%

        steps_matrix = cell2mat(steps_cell)

        steps_matrix = steps_matrix(:,~isnan(sum(steps_matrix,1)))
        
        mu_max = max(steps_matrix,[],2)
        mu_min = min(steps_matrix,[],2)
        sigma_max = max(mu_max-mu_min)
        
        D_guess = std2(steps_matrix)^2/(2*delta_t)
        
        residuals = zeros(max_states,numel(data_cdf));
        
        for num_states = 1:max_states

            beta0 = [ones(1,num_states) / num_states, D_guess*ones(1,num_states)];

            problem = createOptimProblem('lsqcurvefit','x0',beta0,'objective',@(b,x) cdf_modelfun(b,x,delta_t),...
            'lb',zeros(1,num_states*2),'ub',[ones(1,num_states),sigma_max*ones(1,num_states)],'xdata',x_locations,'ydata',data_cdf);

            ms = MultiStart('PlotFcns',@gsplotbestf);
            [beta,errormulti] = run(ms,problem,num_iterations);

            num_D = numel(beta) / 2;

            weights = beta(1:num_D);
            weights = weights / sum(weights);
            Ds = beta(num_D+1:end);

            [Dsort, Is] = sort(Ds,'descend');
            
            D_cell{1,num_states} = Dsort;
            weight_cell{1,num_states} = weights(Is);

            residuals(num_states,:) = data_cdf-cdf_modelfun(beta,x_locations,delta_t);
            
        end
        
        residuals;
        
        ss_residuals = sum(residuals.^2,2)
        
        selected_model = max_states;
        
        for num_states = 1:max_states-1
            k = num_states;
            p = 2;
            n = numel(data_cdf);
            
            F_stat = (ss_residuals(num_states) - ss_residuals(num_states+1)) / p ...
                / (ss_residuals(num_states+1) / ...
                (n-(k+p+1)))
            
            df1 = n - k;
            df2 = n - (k + p);
            
            ss_residuals(num_states) - ss_residuals(num_states+1)
            
%             F_stat = (ss_residuals(num_states) - ss_residuals(num_states+1)) / (df1 - df2) / (ss_residuals(num_states+1) / df2)
            
            P = 1 - fcdf(F_stat,df1,df2)

            if P > F_test_alpha
                selected_model = num_states;
                break
            end
            
        end
        
        selected_model
        
%         return
        
        figure(13423+channel_idx)
        clf
        
            for state_idx = 1:max_states
                
                subplot(2,max_states,state_idx)
                    hold all

                    plot(x_locations,cdf_modelfun([weight_cell{state_idx},D_cell{state_idx}],x_locations,delta_t),'g','LineWidth',5)
                    plot(x_locations,data_cdf,'r-.','LineWidth',3)
                    
                    xlim([min(x_locations), max(x_locations)])

                    set(gca,'xscale','log')
                    
                    xlabel('Displacement (\mum)')
                    ylabel('P(r,\Deltat)')
                    title(['CDF for ' num2str(state_idx) ' states'])

                subplot(2,max_states,max_states+state_idx)
                    hold all
                    
                    residuals = data_cdf-cdf_modelfun([weight_cell{state_idx},D_cell{state_idx}],x_locations,delta_t);
                    
                    xlim([min(x_locations), max(x_locations)])
                    line(xlim,[0 0],'LineStyle','--','LineWidth',3,'Color','k')
                    
                    plot(x_locations,residuals,'r-','LineWidth',3)
                    
                    set(gca,'xscale','log')
                    
                    xlabel('Displacement (\mum)')
                    ylabel('Residuals')    
            end
            
        suptitle(['CDF fitting for: ' options.channel_labels{channel_idx}])
            
        set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')  
        set(gcf,'Color','white')  
            
        figure(13432+channel_idx)
        clf

            for state_idx = 1:max_states
                % Compile conditions

                Ds = D_cell{state_idx}';
                weights = weight_cell{state_idx}';

                subplot(2,max_states,state_idx)

                    bar(Ds)

                    title(['D coeffs for ' num2str(state_idx) ' states'])
                    ylabel('Diff coeff (\mum^2/s)')

                subplot(2,max_states,max_states+state_idx)

                    bar(weights)

                    title(['Weights for ' num2str(state_idx) ' states'])
                    ylabel('Weight (%)')
            end

        suptitle(['Diffusion coefficients for ' num2str(state_idx) ' states'])
            
        set(findall(gcf,'type','text'),'fontSize',14,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',14,'fontWeight','bold')  
        set(gcf,'Color','white')    

    end

    
end

function y = cdf_modelfun(b,x,delta_t)
    %%% Coefficients
    % weights followed by diffusion coefficients

    num_D = numel(b) / 2;
    
    weights = b(1:num_D);
    weights = weights / sum(weights);
    Ds = b(num_D+1:end);
    
    exponentials = zeros(numel(x),num_D);
    
    for D_idx = 1:num_D
        
        exponentials(:,D_idx) = weights(D_idx) .* exp(-x.^2 ./ (4*Ds(D_idx)*delta_t));
    end
    
    y = 1 - sum(exponentials,2);

end

function interaction_results(cell_tracking_results_struct, cell_interaction_results_struct,options)
    num_wells = numel(cell_interaction_results_struct.interactions);
    
    warning('interaction_results: hack where it assumes 2 cell channels')
    
    unique_interaction_counts = zeros(num_wells, 3, 5); % the last is how many histogram bars we'll have
    taulist = cell(1,3);
    
    % wells where there are no <x> cells should not contribute to histogram
    
    for well_idx = 1:num_wells
        
        for UIC_col_idx = 1:numel(options.cell_channels)
            channel_idx = UIC_col_idx + options.cell_channels(1) - 1;
            
            props_struct = cell_tracking_results_struct.linked_object_cells{well_idx,channel_idx};
            
            if numel(props_struct) < 2
                unique_interaction_counts(well_idx,UIC_col_idx,:) = NaN;
                
                if numel(props_struct) < 1
                    unique_interaction_counts(well_idx,3,:) = NaN;
                end
            end
        end
    end
    
    for well_idx = 1:num_wells
        well_interactions = cell_interaction_results_struct.interactions{well_idx};
        
        for interaction_idx = 1:numel(well_interactions)
            int = well_interactions(interaction_idx);
            
            A = int.interacting_flag;
            B = [0 A 0];
            
            % how many unique interactions for these two cells
            changes = B(2:end) - B(1:end-1);
            UIC_count_idx = sum(changes == -1) + 1;
            
            if UIC_count_idx >= size(unique_interaction_counts,3)
                UIC_count_idx = size(unique_interaction_counts,3);
            end
            
            channel_offset = options.cell_channels(1) - 1;
            
            if int.channel1 == int.channel2
                UIC_col = int.channel1 - channel_offset;
            else
                UIC_col = 3;
            end
            
            unique_interaction_counts(well_idx, UIC_col, UIC_count_idx) = unique_interaction_counts(well_idx, UIC_col, UIC_count_idx) + 1;
            
            % how long did each interaction last?
            
            % Find transition points (this is from HMM-Bayes analysis)
            
            transpoints = find(int.interacting_flag(2:end) - int.interacting_flag(1:end-1));
            
            cur_idx = 1;

            for i = 1:numel(transpoints)
                state = int.interacting_flag(transpoints(i));
                
                if state
                    if i == 1
                        taulist{UIC_col} = [taulist{UIC_col}, transpoints(i)-cur_idx+1];
                    else
                        taulist{UIC_col} = [taulist{UIC_col}, transpoints(i)-cur_idx];
                    end
                end

                cur_idx = transpoints(i);
            end

            % Add the last one

            state = int.interacting_flag(end);
            
            if state
                taulist{UIC_col} = [taulist{UIC_col}, numel(int.interacting_flag)-cur_idx];
            end

            if any(taulist{UIC_col} == 0)
                error
            end
        end
    end
    
    UIC_hist_loc = 0:size(unique_interaction_counts,3)-1;
    
    comb_string = [options.channel_labels{options.cell_channels(1)} ' <-> ' options.channel_labels{options.cell_channels(2)} ' interaction'];
    
    UIC_frac_1 = squeeze(unique_interaction_counts(:,1,:)) ./ repmat(sum(squeeze(unique_interaction_counts(:,1,:)),2),[1,size(unique_interaction_counts,3)]);
    UIC_frac_2 = squeeze(unique_interaction_counts(:,2,:)) ./ repmat(sum(squeeze(unique_interaction_counts(:,2,:)),2),[1,size(unique_interaction_counts,3)]);
    UIC_frac_3 = squeeze(unique_interaction_counts(:,3,:)) ./ repmat(sum(squeeze(unique_interaction_counts(:,3,:)),2),[1,size(unique_interaction_counts,3)]);
    
    figure(134287)
    clf
        subplot(1,3,1)
            hold all

%             errorbar(UIC_hist_loc,nanmean(UIC_frac_1,1),nanstd(UIC_frac_1,0,1),'.-','Color','r','LineWidth',3,'MarkerSize',40)
%             errorbar(UIC_hist_loc,nanmean(UIC_frac_2,1),nanstd(UIC_frac_2,0,1),'.-','Color','g','LineWidth',3,'MarkerSize',40)
%             errorbar(UIC_hist_loc,nanmean(UIC_frac_3,1),nanstd(UIC_frac_3,0,1),'.-','Color','b','LineWidth',3,'MarkerSize',40)

            plot(UIC_hist_loc,nanmean(UIC_frac_1,1),'.-','Color','r','LineWidth',3,'MarkerSize',40)
            plot(UIC_hist_loc,nanmean(UIC_frac_2,1),'.-','Color','g','LineWidth',3,'MarkerSize',40)
            plot(UIC_hist_loc,nanmean(UIC_frac_3,1),'.-','Color','b','LineWidth',3,'MarkerSize',40)

            xl = xlim;
            yl = ylim;
            
            xlim([-0.5, 4.5])
            ylim([0, 1])
            
            xlabel('Number of unique interactions')
            ylabel('Fraction of total interactions')
            title('Number of unique cell-cell interactions found per well')
            
            box on
            grid on
            
        subplot(1,3,2)
            hold all

            bar(1, mean(taulist{1}), 'facecolor', 'r'); 
            bar(2, mean(taulist{2}), 'facecolor', 'g'); 
            bar(3, mean(taulist{3}), 'facecolor', 'b');
            
            errorbar(1, mean(taulist{1}), std(taulist{1}), 'k', 'LineWidth',5)%, 'linestyle', 'none');
            errorbar(2, mean(taulist{2}), std(taulist{2}), 'k', 'LineWidth',5)%, 'linestyle', 'none');
            errorbar(3, mean(taulist{3}), std(taulist{3}), 'k', 'LineWidth',5)%, 'linestyle', 'none');
            
%             sigstar(NSPF_groups,NSPF_p_values);
            
            title('Average lifetime of interaction')
            ylabel('Lifetime (frames)')
            xlabel('Cell type / interaction')

            yl = ylim;
            ylim([0, yl(2)])

            set(gca,'XTick',1:3)
            set(gca,'XTickLabel','')
            
        subplot(1,3,3)
            hold all
            
            plot(-1,-1,'.-','Color','r','LineWidth',3,'MarkerSize',40)
            plot(-1,-1,'.-','Color','g','LineWidth',3,'MarkerSize',40)
            plot(-1,-1,'.-','Color','b','LineWidth',3,'MarkerSize',40)
            
            xlim([0 1])
            ylim([0 1])
            axis off
            
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