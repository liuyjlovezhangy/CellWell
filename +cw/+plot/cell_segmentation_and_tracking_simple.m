function movie_file = plot_cell_segmentation_and_tracking(tracking_mode, wells, cell_segmentation_results_struct, cell_tracking_results_struct, is_noise_matrix, make_movies, movie_dir)

    channel_labels = {'K562','CD19 CART','SYTOX'};

    draw_track_len = 15;

    movie_file = [movie_dir '/*.avi'];

    cell_masks_final = cell_segmentation_results_struct.cell_masks_final;
    
    watershed_inputs = cell_segmentation_results_struct.watershed_inputs;
    
    threshold_levels = cell_segmentation_results_struct.threshold_levels;
    thresh_xvals = cell_segmentation_results_struct.thresh_xvals;
    thresh_yvals = cell_segmentation_results_struct.thresh_yvals;

    num_colors = 2^16;
    
    color_space = linspace(0, 1, num_colors)';
    z = zeros(num_colors,1);
    
    r_map = [color_space, z, z];
    g_map = [z, color_space, z];
    b_map = [z, z, color_space];
    
    cmaps = {r_map,g_map,b_map};

    num_signal_channels = size(is_noise_matrix,2);
    
    for well_idx = 1:numel(wells)
        if make_movies
            if tracking_mode
                movie_prefix = 'cell_tracking';
            else
                movie_prefix = 'cell_segmentation';
            end
            writerObj = VideoWriter([movie_dir '/' movie_prefix '_' num2str(well_idx) '_vid.avi']);
            writerObj.FrameRate = 15;
            writerObj.Quality = 100;
            open(writerObj);
        end
        
%         cur_well_img = mat2gray(wells(well_idx).im_well);
        cur_well_img = wells(well_idx).im_well;
        
        if tracking_mode
            figure(132738)
        else
            figure(183299)
        end
        
        for frame_idx = 1:size(cur_well_img,3)
            clf

            for channel_idx = 1:num_signal_channels

                objects = cell_segmentation_results_struct.detected_cell_props_water{well_idx,frame_idx,channel_idx};

                subtightplot(num_signal_channels,4,4*(channel_idx-1) + 1)

                    draw_box(cur_well_img(:,:,frame_idx,channel_idx),0,...
                        objects,well_idx,channel_idx,frame_idx)
                            if channel_idx == 1
                                title('Original image')
                            end
                            
                    ylabel(channel_labels{channel_idx})

                subtightplot(num_signal_channels,4,4*(channel_idx-1) + 2)

                    draw_box(watershed_inputs{well_idx}(:,:,frame_idx,channel_idx),2,...
                        objects,well_idx,channel_idx,frame_idx)
                    if channel_idx == 1
                        title('Watershed input')
                    end

                subtightplot(num_signal_channels,4,4*(channel_idx-1) + 3)

                    draw_box(cell_masks_final{well_idx}(:,:,frame_idx,channel_idx),1,...
                        objects,well_idx,channel_idx,frame_idx)
                    if channel_idx == 1
                        title('Final segmentation')
                    end

                subtightplot(num_signal_channels,4,4*(channel_idx-1) + 4,0.15)

                    hold all

                    xvals = thresh_xvals{well_idx}{frame_idx,channel_idx};
                    yvals = thresh_yvals{well_idx}{frame_idx,channel_idx};
                    
                    xvals(yvals <= 0) = [];
                    yvals(yvals <= 0) = [];
                    
                    plot(xvals,yvals,'-k.','LineWidth',3,'MarkerSize',25)

                    set(gca,'yscale','log')
                    set(gca,'xscale','log')

                    if ~isempty(threshold_levels{well_idx}{frame_idx,channel_idx})
                        line([threshold_levels{well_idx}{frame_idx,channel_idx}.level threshold_levels{well_idx}{frame_idx,channel_idx}.level],ylim,'LineStyle','--','Color','k','LineWidth',3)
                    end

                    xlabel('Normalized pixel intensity')
                    ylabel('Average object area (px)')

                    if channel_idx == 1
                        title('Adaptive thresholding')
                    end

                    pos = get(gca,'OuterPosition');
                    set(gca,'OuterPosition',[pos(1), pos(2) + 0.01, pos(3), pos(4)])

                    box on
                    grid on
                    
                    if ~isempty(xvals)
                        xlim([min(xvals),max(xvals)])
                    end
            end

            suptitle(['Cell segmentation and tracking. Well: ' num2str(well_idx) ' frame: ' num2str(frame_idx)])

            set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
            set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',3)
            set(gcf, 'color', 'white');
            
            set(gcf,'Position',[0 0 1200 1000]);
            fig = gcf;
            fig.PaperPositionMode = 'auto';

            drawnow

            if make_movies
                F=getframe(gcf);

                writeVideo(writerObj,F);
            end
        end

        if make_movies
            close(writerObj);
        end
    end
    
    function draw_box(cur_slice,image_mode,objects,well_idx,channel_idx,frame_idx)
        % image mode
        % 0: image with lines highlighting segmentation
        % 1: segmentation regions with separate colors
        % 2: watershedding input
        
        if isempty(cur_slice)
            return
        end
        
        hold all

        if image_mode == 1
            rgb = label2rgb(cur_slice,'jet',0);
            imagesc(rgb)
        else
            imagesc(cur_slice)
        end

        xlim([1, size(cur_slice,2)])
        ylim([1, size(cur_slice,1)])

        if image_mode == 0
            for obj_idx = 1:numel(objects)

                pts = objects(obj_idx).ConvexHull;
                centroid = objects(obj_idx).Centroid;

                plot(pts(:,1),pts(:,2),'m-','LineWidth',2)
    %             plot(centroid(1),centroid(2),'x','Color','m','MarkerSize',15,'LineWidth',2)

            end
        end
        
        if tracking_mode
            cell_tracks = cell_tracking_results_struct.cell_tracks{well_idx,channel_idx};

            for track_idx = 1:numel(cell_tracks)
                first_draw_track_frame = max([1, frame_idx-draw_track_len]);
                
                track = cell_tracks{track_idx};

%                 track_nan_start_idx = find(~isnan(sum(track,1)),1,'first');
%                 track_nan_end_idx = find(~isnan(sum(track,1)),1,'last');
% 
% %                             padded_track_cmap = [ones(track_nan_start_idx-1,3); track_cmap; ones(size(cur_well_img,3) - track_nan_last_idx,3)];
% 
%                 final_start = max([first_draw_track_frame, track_nan_start_idx]);
%                 final_end = min([frame_idx, track_nan_end_idx]);
% 
%                 track_cmap = spring(numel(final_start:final_end));
                
%                 if isempty(final_start)
%                     continue
%                 end
                
                marker_color = {'y','k','w'};
                line_color = {  'g','m','m';
                                'r','m','m';
                                'g','m','m'};
                
                plot(track(1,frame_idx),track(2,frame_idx),'x','MarkerSize',20,'LineWidth',3,'Color',marker_color{image_mode+1});
                plot(track(1,first_draw_track_frame:frame_idx),track(2,first_draw_track_frame:frame_idx),'-','LineWidth',3,'Color',line_color{channel_idx,image_mode+1})
                
%                 for plot_frame_idx = 2:numel(final_start:final_end)
%                     plot(track(1,plot_frame_idx + final_start - 1),track(2,plot_frame_idx + final_start - 1),'w-','LineWidth',3,'Color',track_cmap(plot_frame_idx,:))
%                 end
                
%                 cline(track(1,final_start:final_end),...
%                     track(2,final_start:final_end),...
%                     zeros(1,final_end-final_start+1),...
%                     1:final_end-final_start,...
%                     track_cmap,...
%                     3);
            end
        end

        axis image
        set(gca,'Ydir','Reverse')
%         axis off
        set(gca,'Color','white')
        set(gca,'XTick',[])
        set(gca,'YTick',[])

        cw.plot.noise_box(is_noise_matrix(well_idx,channel_idx,frame_idx));

        dist = cell_masks_final{well_idx}(:,:,frame_idx,channel_idx);
        dist = dist(:);
        if all(dist==0)
            colormap([0 0 0])
        else
            colormap(cmaps{channel_idx})
        end

        freezeColors
    end
end


