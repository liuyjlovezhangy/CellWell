function test_new_segmentation

    %%% THIS IS FOR 75 um WELLS

%     well_counts = [6,6]; % #rows x #cols
%     min_well_area = 2000;
%     max_well_area = 6000;
%     
%     well_width = 61;
%     well_spacing_width = 14;
%     well_well_dist_width = well_width + well_spacing_width;
%     
%     well_height = 64;
%     well_spacing_height = 11;
%     well_well_dist_height = well_height + well_spacing_height;

%     im = mat2gray(zloadim('C2-CARTonly_75um_BF.tif'));
%     im = mat2gray(zloadim('Timelapse01_BF.tif'));
    
    %%% THIS IS FOR 125 um WELLS

    well_counts = [4,4]; % #rows x #cols
    min_well_area = 8000;
    max_well_area = 11000;
    
    well_width = 97;
    well_spacing_width = 14;
    well_well_dist_width = well_width + well_spacing_width;
    well_array_width = well_width * well_counts(2) + well_spacing_width * (well_counts(2)-1);
    
    well_height = 98;
    well_spacing_height = 13;
    well_well_dist_height = well_height + well_spacing_height;
    well_array_height = well_height * well_counts(1) + well_spacing_height * (well_counts(1)-1);
    
    im = mat2gray(zloadim('CARTonly_125um_BF.tif'));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% init
    
    num_frames = size(im,3);
    
    extrema_extrema_x = repmat([inf,1],[num_frames,1]);
    extrema_extrema_y = repmat([inf,1],[num_frames,1]);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% drift correct
    
%     usfac = 1;
%     
%     multiWaitbar('CloseAll');
%     multiWaitbar('Registering...',0);
%     
%     for frame_idx = 1:num_frames
%         
%         [output, Greg] = dftregistration(fft2(im(:,:,1)),fft2(im(:,:,frame_idx)),usfac);
% 
%         shiftI = round(output(3));
%         shiftJ = round(output(4));
%         
%         im(:,:,frame_idx) = shift_image(im(:,:,frame_idx),shiftI,shiftJ);
%         
%         multiWaitbar('Registering...',frame_idx/num_frames);
%     end
%     
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run segmentation
    
    im_bknd = zeros([size(im(:,:,1)),num_frames]);
    im_noise = zeros([size(im(:,:,1)),num_frames]);
    im_contrast = zeros([size(im(:,:,1)),num_frames]);
    im_edge = zeros([size(im(:,:,1)),num_frames]);
    im_close = zeros([size(im(:,:,1)),num_frames]);
    im_fill = zeros([size(im(:,:,1)),num_frames]);
    im_open = zeros([size(im(:,:,1)),num_frames]);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Segmenting...',0);
    
    for frame_idx = 1:num_frames
        
        cur_frame = im(:,:,frame_idx);
        
        %%% subtract background
        im_bknd(:,:,frame_idx) = cur_frame - imopen(cur_frame,strel('disk',60));
        
        %%% remove noise
        im_noise(:,:,frame_idx) = wiener2(im_bknd(:,:,frame_idx), [5,5]);
        
        %%% contrast adjust
        im_contrast(:,:,frame_idx) = imadjust(im_noise(:,:,frame_idx),[0.1 0.4],[],0.20);
        
        %%% edge detection / segmentation
        im_edge(:,:,frame_idx) = edge(im_contrast(:,:,frame_idx),'canny');

        %%% line dilating / closing      
        sedilate_x = strel('line',3,90);
        sedilate_y = strel('line',3,0);
        seclose = strel('rectangle',[3 3]);
        
        im_close(:,:,frame_idx) = imdilate(im_edge(:,:,frame_idx),sedilate_x);
        im_close(:,:,frame_idx) = imdilate(im_close(:,:,frame_idx),sedilate_y);
        im_close(:,:,frame_idx) = imclearborder(imclose(im_close(:,:,frame_idx),seclose));
        
        %%% fill holes
        im_fill(:,:,frame_idx) = imfill(im_close(:,:,frame_idx),'holes');
        
        %%% opening
        seopen = strel('rectangle',[50,50]);
        im_open(:,:,frame_idx) = imopen(im_fill(:,:,frame_idx),seopen);
        
        multiWaitbar('Segmenting...',frame_idx/num_frames);
    end
    
    im_mask = im_open;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% rotato
    
    % find best angle based on the segmentation mask instead of original,
    % potentially super noisy image
    
    angles = [];
    
    for frame_idx = 1:num_frames
        cur_angle = horizon(im(:,:,frame_idx),0.1,'hough')
        angles = [angles, cur_angle];
    end
    
    angle = mean(angles)
     
    im_rotate = zeros([size(imrotate(im(:,:,1),-angle)), num_frames]);
    im_bknd_rotate = im_rotate;
    im_noise_rotate = im_rotate;
    im_contrast_rotate = im_rotate;
    im_edge_rotate = im_rotate;
    im_close_rotate = im_rotate;
    im_fill_rotate = im_rotate;
    im_open_rotate = im_rotate;
    im_seg_rotate = im_rotate;
    
    multiWaitbar('CloseAll');
    multiWaitbar('Rotating...',0);
    
    for frame_idx = 1:num_frames
        im_rotate(:,:,frame_idx) = imrotate(im(:,:,frame_idx),-angle);
        im_bknd_rotate(:,:,frame_idx) = imrotate(im_bknd(:,:,frame_idx),-angle);
        im_noise_rotate(:,:,frame_idx) = imrotate(im_noise(:,:,frame_idx),-angle);
        im_contrast_rotate(:,:,frame_idx) = imrotate(im_contrast(:,:,frame_idx),-angle);
        im_edge_rotate(:,:,frame_idx) = imrotate(im_edge(:,:,frame_idx),-angle);
        im_close_rotate(:,:,frame_idx) = imrotate(im_close(:,:,frame_idx),-angle);
        im_fill_rotate(:,:,frame_idx) = imrotate(im_fill(:,:,frame_idx),-angle);
        im_open_rotate(:,:,frame_idx) = imrotate(im_open(:,:,frame_idx),-angle);
        
        multiWaitbar('Rotating...',frame_idx/num_frames);
    end
    
    im_mask_rotate = im_open_rotate;
    
    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% infer well centers and generate new masks
    
    multiWaitbar('CloseAll');
    multiWaitbar('Detecting objects...',0);
    
    objects = cell(1,num_frames);
    
    %%% object detection to find extrema of well spacing in this frame
    for frame_idx = 1:num_frames
        im_seg_rotate(:,:,frame_idx) = im_rotate(:,:,frame_idx) .* im_mask_rotate(:,:,frame_idx);
        
        L = bwlabeln(im_seg_rotate(:,:,frame_idx));
        objects{frame_idx} = regionprops(L,im_seg_rotate(:,:,frame_idx),'Extrema','Area','Centroid');
        
        for obj_idx = 1:numel(objects{frame_idx})
            if objects{frame_idx}(obj_idx).Area > min_well_area && objects{frame_idx}(obj_idx).Area < max_well_area
                
                extrema_extrema_x(frame_idx,1) = min([extrema_extrema_x(frame_idx,1), min(objects{frame_idx}(obj_idx).Extrema(:,1))]);
                extrema_extrema_x(frame_idx,2) = max([extrema_extrema_x(frame_idx,2), max(objects{frame_idx}(obj_idx).Extrema(:,1))]);
                
                extrema_extrema_y(frame_idx,1) = min([extrema_extrema_y(frame_idx,1), min(objects{frame_idx}(obj_idx).Extrema(:,2))]);
                extrema_extrema_y(frame_idx,2) = max([extrema_extrema_y(frame_idx,2), max(objects{frame_idx}(obj_idx).Extrema(:,2))]);
            end

        end
        
        multiWaitbar('Detecting objects...',frame_idx/num_frames);
    end
    
    multiWaitbar('CloseAll');
    
    mask_final = zeros(size(im_rotate));
    im_seg_final = mask_final;
    
    %%% repair extrema matrices for frames that didn't segment well
    
    % sanity check
    
    extrema_extrema_x(extrema_extrema_x(:,1) + well_array_width > size(im_rotate,2),:) = NaN;
    extrema_extrema_y(extrema_extrema_y(:,1) + well_array_height > size(im_rotate,1),:) = NaN;
    
    extrema_extrema_x(sum(extrema_extrema_x,2) < 0.9*well_array_width,:) = NaN;
    extrema_extrema_y(sum(extrema_extrema_y,2) < 0.9*well_array_height,:) = NaN;
    
    % outlier detection
    
    alpha = 0.25;
    
    extrema_extrema_x = [deleteoutliers(extrema_extrema_x(:,1),alpha,1), deleteoutliers(extrema_extrema_x(:,2),alpha,1)];
    extrema_extrema_y = [deleteoutliers(extrema_extrema_y(:,1),alpha,1), deleteoutliers(extrema_extrema_y(:,2),alpha,1)];
    
    % sweep extrema list and set NaNs to nearest neighbor extrema
    
    for frame_idx = 2:num_frames
        if isnan(sum(extrema_extrema_x(frame_idx,:))) || isnan(sum(extrema_extrema_y(frame_idx,:)))
            extrema_extrema_x(frame_idx,:) = extrema_extrema_x(frame_idx-1,:);
            extrema_extrema_y(frame_idx,:) = extrema_extrema_y(frame_idx-1,:);
        end
    end
    
    for frame_idx = num_frames-1:-1:1
        if isnan(sum(extrema_extrema_x(frame_idx,:))) || isnan(sum(extrema_extrema_y(frame_idx,:)))
            extrema_extrema_x(frame_idx,:) = extrema_extrema_x(frame_idx+1,:);
            extrema_extrema_y(frame_idx,:) = extrema_extrema_y(frame_idx+1,:);
        end
    end
    
    % make frame-by-frame masks
    
    for frame_idx = 1:num_frames
        
        for row_idx = 1:well_counts(1)
            for col_idx = 1:well_counts(2)
                start_i = round(extrema_extrema_y(frame_idx,1)) + (row_idx-1) * well_well_dist_height;
                end_i = start_i + well_height;
                start_j = round(extrema_extrema_x(frame_idx,1)) + (col_idx-1) * well_well_dist_width;
                end_j = start_j + well_width;

                mask_final(start_i:end_i, start_j:end_j, frame_idx) = 1;
            end 
        end
        
        im_seg_final(:,:,frame_idx) = im_rotate(:,:,frame_idx) .* mask_final(:,:,frame_idx);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plot
        
    figure(134523)
    
    for frame_idx = 1:num_frames
        
        clf
        
        subtightplot(2,5,1)
        hold all
        
            imagesc(im_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Original')
            
        subtightplot(2,5,2)
        hold all
        
            imagesc(im_bknd_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Bknd subtracted')
            
        subtightplot(2,5,3)
        hold all
        
            imagesc(im_noise_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Noise filtered')
            
        subtightplot(2,5,4)
        hold all
        
            imagesc(im_contrast_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Contrast adjusted')
            
        subtightplot(2,5,5)
        hold all
        
            imagesc(im_edge_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Edge detected')
            
        subtightplot(2,5,6)
        hold all
        
            imagesc(im_close_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Closed, border cleared')
            
        subtightplot(2,5,7)
        hold all
        
            imagesc(im_fill_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Filled holes')
            
        subtightplot(2,5,8)
        hold all
        
            imagesc(im_open_rotate(:,:,frame_idx))

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('Opened')
        
        subtightplot(2,5,9)
        hold all
        
            imagesc(im_seg_rotate(:,:,frame_idx))
            
            for obj_idx = 1:numel(objects{frame_idx})

                extrema = objects{frame_idx}(obj_idx).Extrema;
                
                xmin = min(extrema(:,1));
                xmax = max(extrema(:,1));
                ymin = min(extrema(:,2));
                ymax = max(extrema(:,2));
                
                pts = [ xmin ymin;...
                        xmin ymax;...
                        xmax ymax;...
                        xmax ymin;...
                        xmin ymin;];
                
                centroid = objects{frame_idx}(obj_idx).Centroid;
                
                if objects{frame_idx}(obj_idx).Area < min_well_area || objects{frame_idx}(obj_idx).Area > max_well_area
                    plot(pts(:,1),pts(:,2),'r','LineWidth',2)
                else
                    plot(pts(:,1),pts(:,2),'g','LineWidth',2)
                    plot(centroid(1),centroid(2),'+','Color','r','MarkerSize',10,'LineWidth',3)
                end

            end

            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,1)-5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,2)+5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,1)-5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,2)+5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Classified')
            
        subtightplot(2,5,10)
        hold all
        
            imagesc(im_seg_final(:,:,frame_idx))
            
            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,1)-5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,2)+5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,1)-5,extrema_extrema_x(frame_idx,1)-5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');
            line([extrema_extrema_x(frame_idx,2)+5,extrema_extrema_x(frame_idx,2)+5],[extrema_extrema_y(frame_idx,1)-5,extrema_extrema_y(frame_idx,2)+5],'LineStyle','-','LineWidth',3,'Color','b');

            colormap gray
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')

            title('Inferred')
            
        suptitle(['Frame ' num2str(frame_idx) ' of ' num2str(num_frames)])
            
        set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');
        drawnow
    end
end
