function cell_segmentation_and_tracking_circle(tracking_mode, wells, cell_segmentation_results_struct, cell_tracking_results_struct, is_noise_matrix, options, make_movies, movie_dir)

%%%%%% NOTE: THIS FILE ADJUSTS GAMMA ON YOUR MOVIES IN ORDER TO BE ABLE TO
%%%%%% PROPERLY SEE DARK OBJECTS

    plopts = options.plot_options;

    warning('cell_segmentation_and_tracking_circle assumes that the cell channels are right next to each other')

    num_wells = numel(wells);
    num_frames = size(wells(1).im_well,3);
    num_channels = size(wells(1).im_well,4);
    
    channel_labels = options.channel_labels;

    cell_channels = sort(options.cell_channels);
    
    draw_track_len = 15;
    gamma_adjust = 0.3;

    cell_masks_final = cell_segmentation_results_struct.cell_masks;
    
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
    
    if tracking_mode
        movie_prefix = 'cell_tracking';
    else
        movie_prefix = 'cell_segmentation';
    end
    
    delete([movie_dir '/' movie_prefix '_*_vid.avi']);
    
    for well_idx = 1:num_wells
        
        cell_counts = zeros(1,num_channels);
        
        for frame_idx = 1:num_frames
            for channel_idx = cell_channels
                objects = cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,channel_idx};

                cell_counts(channel_idx) = max([cell_counts(channel_idx), numel(objects)]);
            end
        end
        
        if any(cell_counts(cell_channels) < plopts.min_cells)
            continue
        end
        
        if make_movies
            writerObj = VideoWriter([movie_dir '/' movie_prefix '_' num2str(well_idx) '_vid.avi']);
            writerObj.FrameRate = 10;
            writerObj.Quality = 100;
            open(writerObj);
        end
        
        cur_well_img = wells(well_idx).im_well;
        
        if tracking_mode
            figure(132738)
        else
            figure(183299)
        end
        
        for frame_idx = 1:num_frames
            clf
            
            for channel_idx = cell_channels

                objects = cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,channel_idx};
    
                subtightplot(numel(cell_channels),2,2*(channel_idx-cell_channels(1)) + 1)
                    hold all

                    draw_box(imadjust(mat2gray(cur_well_img(:,:,frame_idx,channel_idx)),[],[],gamma_adjust),0,...
                        objects,well_idx,channel_idx,frame_idx)
                        if channel_idx == 1
                            title('Original image')
                        end
                        
                    ylabel(channel_labels{channel_idx})
                    
                subtightplot(numel(cell_channels),2,2*(channel_idx-cell_channels(1)) + 2)

                    draw_box(cell_masks_final{well_idx}(:,:,frame_idx,channel_idx),1,...
                        objects,well_idx,channel_idx,frame_idx)
                    if channel_idx == 1
                        title('Final segmentation')
                    end
            end

            suptitle(['Cell segmentation and tracking. Well: ' num2str(well_idx) ' frame: ' num2str(frame_idx)])

            set(findall(gcf,'type','text'),'fontSize',16,'fontWeight','bold')
            set(findall(gcf,'type','axes'),'fontSize',16,'fontWeight','bold','LineWidth',5)
            set(gcf, 'color', 'white');
            
%             set(gcf,'Position',[0 0 1200 1000]);
%             fig = gcf;
%             fig.PaperPositionMode = 'auto';

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
            
            rgb = label2rgb(cur_slice,'jet',0.8*[1 1 1]);
            imagesc(rgb)
            
            box on
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
                
                marker_color = {'y','k','w'};
                line_color = {  'g','m','m';
                                'r','m','m';
                                'g','m','m'};
                
                plot(track(1,frame_idx),track(2,frame_idx),'x','MarkerSize',10,'LineWidth',3,'Color',marker_color{image_mode+1});
                plot(track(1,first_draw_track_frame:frame_idx),track(2,first_draw_track_frame:frame_idx),'-','LineWidth',3,'Color',line_color{channel_idx,image_mode+1})

            end
        end

        axis image
        set(gca,'Ydir','Reverse')
%         axis off
        
        set(gca,'Color','white')
        set(gca,'XTick',[])
        set(gca,'YTick',[])

%         cw.plot.noise_box(is_noise_matrix(well_idx,channel_idx,frame_idx));

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


