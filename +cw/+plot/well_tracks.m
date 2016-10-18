function movie_file = plot_well_tracks(im, im_shifted, seg_results_struct, well_tracking_results_struct, make_movies, movie_dir )

    movie_file = [movie_dir '/well_tracking.avi'];

    if make_movies
        writerObj = VideoWriter(movie_file);
        writerObj.FrameRate = 18;
        writerObj.Quality = 100;
        open(writerObj);
    end

    figure(78193)
    clf

    for frame_idx = 1:size(im,3)
        bf_seg = im(:,:,frame_idx,end) .* seg_results_struct.frame(frame_idx).im_mask_final_corrected;
        
        bf_seg_shifted = im_shifted(:,:,frame_idx,end) .* well_tracking_results_struct.mask;

        subplot(1,2,1)
            hold all

            imagesc(bf_seg)
            colormap gray

            for track_idx = 1:numel(well_tracking_results_struct.well_tracks)
                track = well_tracking_results_struct.well_tracks{track_idx};

                plot(track(1,1:frame_idx),track(2,1:frame_idx),'g-','LineWidth',3)
                plot(track(1,frame_idx),track(2,frame_idx),'or','LineWidth',2)
            end

            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Pre-drift correction')

        subplot(1,2,2)
            hold all

            imagesc(bf_seg_shifted)
            colormap gray

            for track_idx = 1:numel(well_tracking_results_struct.well_tracks_shifted)
                track = well_tracking_results_struct.well_tracks_shifted{track_idx};

                plot(track(1,1:frame_idx),track(2,1:frame_idx),'g-','LineWidth',3)
                plot(track(1,frame_idx),track(2,frame_idx),'or','LineWidth',2)
            end

            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Post-drift correction with well size correction')

        suptitle(['Well tracking and drift correction; frame ' num2str(frame_idx) ' of ' num2str(size(im,3))])

        set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');

    %     set(gcf,'Position',[0 0 1400 700]);
    %     fig = gcf;
    %     fig.PaperPositionMode = 'auto';

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

