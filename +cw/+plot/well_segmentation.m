function movie_file = plot_well_segmentation( seg_results_struct, im, otsu_idcs, make_movies, movie_dir )

    movie_file = [movie_dir '/segmentation.avi'];

    if make_movies
        writerObj = VideoWriter(movie_file);
        writerObj.FrameRate = 5;
        writerObj.Quality = 100;
        open(writerObj);
    end
    
    for frame_idx = 1:size(im,3)
        
        im_bf_slice = im(:,:,frame_idx,end);
        
        bf_seg = im_bf_slice .* seg_results_struct.frame(frame_idx).im_mask_final_corrected;
        
        %%% Plotting

        figure(81322)
        clf

        subplot(2,3,1)
            hold all

            imagesc(im_bf_slice)

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Original image')

        subplot(2,3,2)
            hold all
            hist(im_bf_slice(:),2000)

            xlabel('Normalized pixel intensity')
            ylabel('# of pixels')
            title('Otsu thresholding')

            for idcs_idx = 1:numel(otsu_idcs)
                line([otsu_idcs(idcs_idx) otsu_idcs(idcs_idx)],ylim,'LineStyle','--','Color','r','LineWidth',3)
            end

            axis tight

        subplot(2,3,3)
            hold all

            imagesc(seg_results_struct.frame(frame_idx).im_mask_thresh)

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Otsu thresholded')

        subplot(2,3,4)
            hold all

            imagesc(seg_results_struct.frame(frame_idx).mask)

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Final mask')

        subplot(2,3,5)
            hold all

            imagesc(im_bf_slice)

            objects = seg_results_struct.frame(frame_idx).all_objects;
            
            for obj_idx = 1:numel(objects)

                extrema = objects(obj_idx).Extrema;
                
                xmin = min(extrema(:,1)) - seg_results_struct.extra_border_x;
                xmax = max(extrema(:,1)) + seg_results_struct.extra_border_x;
                ymin = min(extrema(:,2)) - seg_results_struct.extra_border_y;
                ymax = max(extrema(:,2)) + seg_results_struct.extra_border_y;
                
                pts = [ xmin ymin;...
                        xmin ymax;...
                        xmax ymax;...
                        xmax ymin;...
                        xmin ymin;];
                
                centroid = objects(obj_idx).Centroid;

                if ~well_criterion(objects(obj_idx))
                    plot(pts(:,1),pts(:,2),'r','LineWidth',2)
                else
                    plot(pts(:,1),pts(:,2),'g','LineWidth',2)
                    plot(centroid(1),centroid(2),'+','Color','r','MarkerSize',10,'LineWidth',3)
                end

            end

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Detected wells')

        subplot(2,3,6)
            hold all

            imagesc(bf_seg)
            
            objects = seg_results_struct.frame(frame_idx).all_objects;

            for obj_idx = 1:numel(objects)

                extrema = objects(obj_idx).Extrema;
                
                xmin = min(extrema(:,1)) - seg_results_struct.extra_border_x;
                xmax = max(extrema(:,1)) + seg_results_struct.extra_border_x;
                ymin = min(extrema(:,2)) - seg_results_struct.extra_border_y;
                ymax = max(extrema(:,2)) + seg_results_struct.extra_border_y;
                
                pts = [ xmin ymin;...
                        xmin ymax;...
                        xmax ymax;...
                        xmax ymin;...
                        xmin ymin;];
                    
                centroid = objects(obj_idx).Centroid;

                if ~well_criterion(objects(obj_idx))
                    plot(pts(:,1),pts(:,2),'r','LineWidth',2)
                else
                    plot(pts(:,1),pts(:,2),'g','LineWidth',2)
                    plot(centroid(1),centroid(2),'+','Color','r','MarkerSize',10,'LineWidth',3)
                end

            end

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Segmented wells')

        suptitle(['Well detection performance; Frame: ' num2str(frame_idx) ' out of ' num2str(size(im,3))])

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

