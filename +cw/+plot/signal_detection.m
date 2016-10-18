function movie_file = plot_signal_detection( well_data,detection_images,is_noise_matrix,noise_threshold_mean,noise_threshold_stdev,make_movies,movie_dir )

    movie_file = [movie_dir '/*.avi'];

    well_box_colors = [0 1 0; 1 0 0];
    
    num_colors = 2^16;
    
    color_space = linspace(0, 1, num_colors)';
    z = zeros(num_colors,1);
    
    r_map = [color_space, z, z];
    g_map = [z, color_space, z];
    b_map = [z, z, color_space];
    
    cmaps = {r_map,g_map,b_map};
    
    for well_idx = 1:numel(well_data)
        
        if make_movies
            writerObj = VideoWriter([movie_dir '/well_by_well_signal_detection_' num2str(well_idx) '_vid.avi']);
            writerObj.FrameRate = 20;
            writerObj.Quality = 100;
            open(writerObj);
        end

        cur_well_im = well_data(well_idx).im_well;
        cur_detection_im = detection_images{well_idx};
        
        for frame_idx = 1:size(cur_well_im,3)
            
            figure(10389)
            clf
        
            for channel_idx = 1:size(cur_detection_im,4)
            
                frame_filter_dist = cur_detection_im(:,:,frame_idx,channel_idx);
                frame_filter_dist = frame_filter_dist(:)';
                frame_filter_dist = frame_filter_dist(frame_filter_dist ~= 0);

                subtightplot(size(cur_detection_im,4),3,3*(channel_idx-1) + 1,0.01)
                    hold all

                    imagesc(cur_well_im(:,:,frame_idx,channel_idx))

                    axis image
                    axis tight
                    set(gca,'Ydir','Reverse')
                    axis off

                    plot_noise_box(is_noise_matrix(well_idx,channel_idx,frame_idx));

                    if channel_idx==1
                        title('Raw image')
                    end
                    
                    ylabel(['Channel: ' num2str(channel_idx)])
                    colormap(cmaps{channel_idx})
                    
                    freezeColors

                subtightplot(size(cur_detection_im,4),3,3*(channel_idx-1) + 2,0.01)

                    hold all

                    imagesc(cur_detection_im(:,:,frame_idx,channel_idx))

                    axis image
                    axis tight
                    set(gca,'Ydir','Reverse')
                    axis off

                    plot_noise_box(is_noise_matrix(well_idx,channel_idx,frame_idx));

                    if channel_idx==1
                        title('Filtered image')
                    end
                    
                    dist = cur_detection_im(:,:,frame_idx,channel_idx);
                    dist = dist(:);
                    if all(dist==0)
                        colormap([0 0 0])
                    else
                        colormap(cmaps{channel_idx})
                    end
                    
                    freezeColors

                subtightplot(size(cur_detection_im,4),3,3*(channel_idx-1) + 3,0.1)

                    hold all

                    [frame_filter_hist,x0] = hist(frame_filter_dist,15);
                    plot(x0,frame_filter_hist,'Color','k','LineStyle','-','LineWidth',3)

                    xlim([0, max([frame_filter_dist, noise_threshold_mean])])
                    axis tight

                    line([mean(frame_filter_dist) mean(frame_filter_dist)],ylim,'LineStyle','--','Color','m','LineWidth',3)
                    line([noise_threshold_mean noise_threshold_mean],ylim,'LineStyle','-','Color',well_box_colors(is_noise_matrix(well_idx,channel_idx,frame_idx)+1,:),'LineWidth',3)

                    if channel_idx==1
                        title('Distribution of filter values')
                    end
                    xlabel('Filter value (AU)')
                    ylabel('# of pixels')
                    box on
                    grid on
                    
                    text(1,1,['\mu: ' sprintf('%5.2f',mean(frame_filter_dist)) '/' num2str(noise_threshold_mean) ' || \sigma: ' sprintf('%5.2f',std(frame_filter_dist)) '/' num2str(noise_threshold_stdev)],...
                        'Color',well_box_colors(is_noise_matrix(well_idx,channel_idx,frame_idx)+1,:),...
                        'Units','Normalized','BackgroundColor','w','HorizontalAlignment','right')
                    
%                     is_noise_matrix(well_idx,channel_id
            end
            
            suptitle(['Filtering wells with no cell signal. Well: ' num2str(well_idx) ' frame: ' num2str(frame_idx)])

            set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
            set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
            set(gcf, 'color', 'white');
            
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

end

