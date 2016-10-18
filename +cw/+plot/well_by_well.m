function movie_file = plot_well_by_well( well_tracking_results_struct,make_movies,movie_dir )

    movie_file = [movie_dir '/*.avi'];

    num_colors = 2^16;
    
    color_space = linspace(0, 1, num_colors)';
    z = zeros(num_colors,1);
    
    r_map = [color_space, z, z];
    g_map = [z, color_space, z];
    b_map = [z, z, color_space];

    for well_idx = 1:numel(well_tracking_results_struct.wells)
        if make_movies
            writerObj = VideoWriter([movie_dir '/well_by_well_' num2str(well_idx) '_vid.avi']);
            writerObj.FrameRate = 12;
            writerObj.Quality = 100;
            open(writerObj);
        end
        
        im_well = well_tracking_results_struct.wells(well_idx).im_well;

        % Fluorescence image for a particular channel for this well, normed
        % to the max intensity encountered for that channel in the well
        
        channel_1 = im_well(:,:,:,1);
        channel_2 = im_well(:,:,:,2);
        channel_3 = im_well(:,:,:,3);
        channel_4 = im_well(:,:,:,4);

        channel_1_norm = mat2gray(channel_1);
        channel_2_norm = mat2gray(channel_2);
        channel_3_norm = mat2gray(channel_3);

        % The total fluorescent signal in a frame, normalized to max signal

        integrated_channel_1 = squeeze(sum(sum(channel_1,1),2));
        integrated_channel_2 = squeeze(sum(sum(channel_2,1),2));
        integrated_channel_3 = squeeze(sum(sum(channel_3,1),2));

        integrated_channel_1_norm = mat2gray(integrated_channel_1);
        integrated_channel_2_norm = mat2gray(integrated_channel_2);
        integrated_channel_3_norm = mat2gray(integrated_channel_3);

        % draw it

        for frame_idx = 1:size(im_well,3)
            
            channel_1_cur = mat2gray(channel_1(:,:,frame_idx));
            channel_2_cur = mat2gray(channel_2(:,:,frame_idx));
            channel_3_cur = mat2gray(channel_3(:,:,frame_idx));
            channel_4_cur = mat2gray(channel_4(:,:,frame_idx));
            
            figure(83823)
            clf

                % plot the fluorescence images

                subplot(2,4,1)
                    hold all

                    imagesc(channel_1_cur)

                    colormap(g_map)
                    axis image
                    axis off
                    set(gca,'Ydir','Reverse')

                    title('Cell channel 1')

                    freezeColors

                subplot(2,4,2)
                    hold all

                    imagesc(channel_2_cur)

                    colormap(b_map)
                    axis image
                    axis off
                    set(gca,'Ydir','Reverse')

                    title('Cell channel 2')

                    freezeColors

                subplot(2,4,3)
                    hold all

                    imagesc(channel_3_cur)

                    colormap(r_map)
                    axis image
                    axis off
                    set(gca,'Ydir','Reverse')

                    title('Readout channel')

                    freezeColors

                subplot(2,4,4)
                    hold all

                    imagesc(channel_4_cur)

                    colormap gray
                    axis image
                    axis off
                    set(gca,'Ydir','Reverse')

                    title('Brightfield Channel')

                    freezeColors

                % Fluorescence intensity integration over time

                subplot(2,4,5)
                    hold all

                    plot(1:frame_idx,integrated_channel_1_norm(1:frame_idx),'g-','LineWidth',3)

                    xlabel('Time [frame]')
                    ylabel('Integr. intensity')

                    xlim([0,size(im_well,3)])
                    ylim([min(integrated_channel_1_norm), max(integrated_channel_1_norm)])
                    
                    title('Cell channel 1')

                subplot(2,4,6)
                    hold all

                    plot(1:frame_idx,integrated_channel_2_norm(1:frame_idx),'b-','LineWidth',3)

                    xlabel('Time [frame]')
                    ylabel('Integr. intensity')

                    xlim([0,size(im_well,3)])
                    ylim([min(integrated_channel_2_norm), max(integrated_channel_2_norm)])
                    
                    title('Cell channel 2')

                subplot(2,4,7)
                    hold all

                    plot(1:frame_idx,integrated_channel_3_norm(1:frame_idx),'r-','LineWidth',3)

                    xlabel('Time [frame]')
                    ylabel('Integr. intensity')

                    xlim([0,size(im_well,3)])
                    ylim([min(integrated_channel_3_norm), max(integrated_channel_3_norm)])
                    
                    title('Readout channel')

            % plot location

                subplot(2,4,8)
                    hold all

                    imagesc(well_tracking_results_struct.im_shifted(:,:,frame_idx,end))

                    well_xmin = well_tracking_results_struct.wells(well_idx).left_boundary;
                    well_xmax = well_tracking_results_struct.wells(well_idx).right_boundary;
                    well_ymin = well_tracking_results_struct.wells(well_idx).bottom_boundary;
                    well_ymax = well_tracking_results_struct.wells(well_idx).top_boundary;
                    
                    pts = [ well_xmin well_ymin;...
                            well_xmin well_ymax;...
                            well_xmax well_ymax;...
                            well_xmax well_ymin;...
                            well_xmin well_ymin;];
                    
                    centroid = well_tracking_results_struct.wells(well_idx).track_shifted(:,frame_idx);

                    plot(pts(:,1),pts(:,2),'g','LineWidth',2)
                    plot(centroid(1),centroid(2),'+','Color','r','MarkerSize',10,'LineWidth',3)

                    colormap gray
                    axis image
                    axis off
                    set(gca,'Ydir','Reverse')

                    title('Well location')

            suptitle(['Cropped well image over time; well ' num2str(well_idx) ' of ' num2str(numel(well_tracking_results_struct.wells)) '; frame ' num2str(frame_idx) ' of ' num2str(size(im_well,3))])

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



