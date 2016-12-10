function cell_overlays(wells, signal_detection_results_struct, cell_segmentation_results_struct, cell_tracking_results_struct, options, make_movies, movie_dir)

%     warning('cell_overlays currently only works for a 4 channel movie with signal / cell channels / brightfield')

    delete([movie_dir '/cell_overlays_*_vid.avi'])

    is_noise_matrix = signal_detection_results_struct.is_noise_matrix;
    
    draw_track_len = 15;

    cell_masks_final = cell_segmentation_results_struct.cell_masks;

    num_colors = 2^16;
    
    color_space = linspace(0, 1, num_colors)';
    z = zeros(num_colors,1);
    
    r_map = [color_space, z, z];
    g_map = [z, color_space, z];
    b_map = [z, z, color_space];
    
    cmaps = {r_map,g_map,b_map};
    
    color_strings = {'r','g','b'};
    
    num_wells = numel(wells);
    num_frames = size(wells(1).im_well,3);
    
    figure(13425)
    clf
    
        hold all
        
        plot(-1,-1,'-','Color','r','LineWidth',8)
        plot(-1,-1,'-','Color','g','LineWidth',8)
%         plot(-1,-1,'-','Color','b','LineWidth',8)
        
        xlim([0 1])
        ylim([0 1])
        axis off
        
        legend(options.channel_labels(options.cell_channels))
        
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');

    drawnow
    
    for well_idx = 1:num_wells
        if make_movies
            writerObj = VideoWriter([movie_dir '/cell_overlays_' num2str(well_idx) '_vid.avi']);
            writerObj.FrameRate = 8;
            writerObj.Quality = 100;
            open(writerObj);
        end
        
        cur_well_im = wells(well_idx).im_well;
        
        for frame_idx = 1:num_frames;
            objects = [];
            for channel_idx = options.cell_channels
                objects = [objects, cell_segmentation_results_struct.detected_cell_props{well_idx,frame_idx,channel_idx}];
            end
            
            figure(1833)
            clf

                well_scaled = squeeze(cur_well_im(:,:,frame_idx,:));
                % change channel brightnesses
                
%                 for channel_idx = 1:size(well_scaled,3)
%                     well_scaled(:,:,channel_idx) = mat2gray(well_scaled(:,:,channel_idx));
%                 end
                
                well_scaled(:,:,options.signal_channels) = 0;

                subtightplot(1,2,1)
                    hold all
                    draw_box(well_scaled,0,objects,well_idx,frame_idx)
                    
                subtightplot(1,2,2)
                    hold all
                    draw_box(squeeze(cell_masks_final{well_idx}(:,:,frame_idx,:)),1,objects,well_idx,frame_idx)
                    
            suptitle(['Multi-channel overlay. Well: ' num2str(well_idx) ' frame: ' num2str(frame_idx)])

            set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
            set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
            set(gcf, 'color', 'white');

            set(gcf,'Position',[0 0 1400 700]);
            fig = gcf;
%             fig.PaperPositionMode = 'auto';
            
            drawnow
            
%             pause

            if make_movies
                F=getframe(gcf);

                writeVideo(writerObj,F);
            end
        end
    end

    if make_movies
        close(writerObj);
    end
    
    function draw_box(cur_slice,image_mode,objects,well_idx,frame_idx)
        % image mode
        % 0: image with lines highlighting segmentation
        % 1: segmentation regions with separate colors
        
        if isempty(cur_slice)
            return
        end
        
        hold all

        xlim([1, size(cur_slice,2)])
        ylim([1, size(cur_slice,1)])
        
        if image_mode == 1
%             cmap = [0 0 0;1,0,0;0,1,0;0,0,1;1,0,1];
            cmap = [0 0 0;1,0,0;0,1,0;0,0,1];
            
            combined_layers = zeros(size(cur_slice,1),size(cur_slice,2),size(cur_slice,3));
            combined_layers(:,:,1) = 1;
            
            border_layer = zeros(size(cur_slice,1),size(cur_slice,2));
            objects_layer = zeros(size(cur_slice,1),size(cur_slice,2),size(cur_slice,3));
            
            for channel_box_idx = options.cell_channels
                layer = cur_slice(:,:,channel_box_idx);          
                
                border_layer(layer == 1) = 1;
                objects_layer(:,:,channel_box_idx) = layer > 1;
                
            end
            
            for channel_box_idx = options.cell_channels
                combined_layers(:,:,channel_box_idx) = (channel_box_idx - options.cell_channels(1) + 1)*objects_layer(:,:,channel_box_idx);
            end
            
            combined_layers = sum(combined_layers,3);

            combined_layers(logical(border_layer)) = 0;
            
            rgb = label2rgb(combined_layers,cmap, [.7 .7 .7]);
            imshow(rgb)

        else
            img = cur_slice(:,:,[options.cell_channels,options.signal_channels]);%cat(3,cur_slice(:,:,[options.cell_channels]),zeros(size(cur_slice,1),size(cur_slice,2)));
            
            imagesc(img)
        end

        if image_mode == 0
            for obj_idx = 1:numel(objects)
                pts = objects(obj_idx).ConvexHull;

                plot(pts(:,1),pts(:,2),'m-','LineWidth',2)

            end
        end
        
        for channel_box_idx = options.cell_channels
            cell_tracks = cell_tracking_results_struct.cell_tracks{well_idx,channel_box_idx};

            for track_idx = 1:numel(cell_tracks)
                first_draw_track_frame = max([1, frame_idx-draw_track_len]);
                
                track = cell_tracks{track_idx};
                
                marker_color = {'y','k','w'};
                line_color = {  'w','w';
                                'w','w';
                                'w','w';
                                };
                            
                if frame_idx <= size(track,2)
                    plot(track(1,frame_idx),track(2,frame_idx),'x','MarkerSize',20,'LineWidth',3,'Color',marker_color{image_mode+1});
                    plot(track(1,first_draw_track_frame:frame_idx),track(2,first_draw_track_frame:frame_idx),'-','LineWidth',3,'Color',line_color{channel_box_idx,image_mode+1})
                end
            end
        end

        axis image
        set(gca,'Ydir','Reverse')
        axis off

%         cw.plot.noise_box(~any(squeeze(~is_noise_matrix(well_idx,:,frame_idx))));

        freezeColors
    end
end
    
    

