function [well_segmentation_results_struct, im_seg_final] = segment_wells_otsu2(im, options)

    propts = options.processing_options;

    im_bf = im(:,:,:,options.bf_channel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    well_well_dist_width = options.well_width + options.well_spacing_width;
    well_array_width = options.well_width * options.well_counts(2) + options.well_spacing_width * (options.well_counts(2)-1);
    
    well_well_dist_height = options.well_height + options.well_spacing_height;
    well_array_height = options.well_height * options.well_counts(1) + options.well_spacing_height * (options.well_counts(1)-1);

    min_well_area = 0.8 * options.well_width * options.well_height;
    max_well_area = 1.2 * options.well_width * options.well_height;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% init
    
    num_frames = size(im,3);
    num_channels = size(im,4);
    
    extrema_extrema_x = repmat([inf,1],[num_frames,1]);
    extrema_extrema_y = repmat([inf,1],[num_frames,1]);
    
    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run segmentation
    
    im_bknd = zeros(size(im_bf));
    im_noise = im_bknd;
    im_contrast = im_bknd;
    im_edge = im_bknd;
    im_close = im_bknd;
    im_fill = im_bknd;
    im_open = im_bknd;
    im_seg = im_bknd;
    
    multiWaitbar('CloseAll');
    multiWaitbar('Segmenting...',0);
    
    for frame_idx = 1:num_frames
        
        cur_frame = im_bf(:,:,frame_idx);
        
        %%% uneven illumination && background subtract
        im_bknd(:,:,frame_idx) = imtophat(cur_frame,strel('disk',20));
%         im_bknd(:,:,frame_idx) = cur_frame;
        %%% subtract background
%         im_bknd(:,:,frame_idx) = cur_frame - imopen(cur_frame,strel('disk',60));
        
        %%% remove noise
        im_noise(:,:,frame_idx) = wiener2(im_bknd(:,:,frame_idx), [5,5]);
        
        %%% contrast adjust

        [r, c]=size(im_noise(:,:,frame_idx));

%         im_contrast(:,:,frame_idx) = mat2gray(homofil(im_noise(:,:,frame_idx),d,r,c,order));
        
        im_contrast(:,:,frame_idx) = imadjust(im_noise(:,:,frame_idx),[],[],0.5);
        
        levels_prev = [];
        
        [level, level_low, level_high, ...
            cur_thresh_xvals, cur_thresh_yvals] = ...
            cw.process.adaptiveThresholdFinder(cur_frame, propts.cseg_adaptive_thresh_scale,...
            propts.cseg_tracking_params, levels_prev);
        
        im_thresh = threshold3D(im_contrast(:,:,frame_idx), level);
        
        
        figure(123523)
            clf
                subplot(1,2,1)
                    hold all

                    xvals = cur_thresh_xvals;
                    yvals = cur_thresh_yvals;

                    xvals(yvals <= 0) = [];
                    yvals(yvals <= 0) = [];

                    plot(xvals,yvals,'-k.','LineWidth',3,'MarkerSize',25)

                    set(gca,'yscale','log')
                    set(gca,'xscale','log')

                    if ~isempty(level)
                        line([level level],ylim,'LineStyle','--','Color','k','LineWidth',3)
                    end

                    xlabel('Normalized pixel intensity')
                    ylabel('Average object area (px)')

                        title('Adaptive thresholding')

                    pos = get(gca,'OuterPosition');
                    set(gca,'OuterPosition',[pos(1), pos(2) + 0.01, pos(3), pos(4)])

                    box on
                    grid on

                    if ~isempty(xvals)
                        xlim([min(xvals),max(xvals)])
                    end

                subplot(1,2,2)
                    hold all
                    imagesc(im_thresh)

                    axis image
                    set(gca,'Ydir','Reverse')
                    axis off 

                    colormap gray
        
        error
        
        %%% edge detection / segmentation
%         im_edge(:,:,frame_idx) = edge(im_contrast(:,:,frame_idx),'canny');
%         x = 1;
%         mu = -1;
%         sigma = 20;
%         sigma_width = 10;

        otsu_max = 2;
        num_otsu = 5;
        
%         data_out = plus_filt2D(im_noise(:,:,frame_idx),x,mu,sigma,sigma_width);
% 
%         data_out = mat2gray(data_out);
%         
%         [otsu_image,~,thresh] = otsu(data_out,num_otsu);
        otsu_input = im_contrast(:,:,frame_idx);
        otsu_input(otsu_input==0) = NaN;
        
        [otsu_out,~,thresh] = otsu(otsu_input,num_otsu);
        
        edge_input = im_contrast(:,:,frame_idx);
        edge_input(otsu_out(:) ~= otsu_max) = 0;
        
%         
%         im_edge(:,:,frame_idx) = data_out < thresh(1);

        im_edge(:,:,frame_idx) = edge_input;%edge(edge_input,'canny');

        %%% line dilating / closing      
        sedilate_x = strel('line',3,90);
        sedilate_y = strel('line',3,0);
        seclose = strel('rectangle',[6 6]);
        
        im_close(:,:,frame_idx) = imdilate(im_edge(:,:,frame_idx),sedilate_x);
        im_close(:,:,frame_idx) = imdilate(im_close(:,:,frame_idx),sedilate_y);
        im_close(:,:,frame_idx) = imclose(im_close(:,:,frame_idx),seclose);
        
        %%% fill holes
        im_fill(:,:,frame_idx) = imfill(imclearborder(im_close(:,:,frame_idx)),'holes');
        
        %%% opening
        seopen = strel('rectangle',[50,50]);
        im_open(:,:,frame_idx) = imopen(im_fill(:,:,frame_idx),seopen);
        
%         im_open(:,:,frame_idx) = im_fill(:,:,frame_idx);

        if options.processing_options.wseg_debug
            figure(14332)
            clf
            
                subtightplot(2,4,1)
                    hold all
                    
                    imagesc(im_bknd(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(2,4,2)
                    hold all
                    
                    imagesc(im_noise(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(2,4,3)
                    hold all
                    
                    imagesc(im_contrast(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(2,4,4)
                    hold all
                    a = im_contrast(:,:,frame_idx);
                    hist(a(:),100)
                    
                    for thresh_idx = 1:numel(thresh)
                        line([thresh(thresh_idx) thresh(thresh_idx)], ylim,'Color','k','LineWidth',3,'LineStyle','--')
                    end
                    
%                     imagesc(otsu_image)
%                     
%                     colormap gray
%                     
%                     axis image
%                     set(gca,'Ydir','Reverse')
% 
%                     set(gca,'Color','white')
%                     set(gca,'XTick',[])
%                     set(gca,'YTick',[])
                    
                subtightplot(2,4,5)
                    hold all
                    
                    imagesc(otsu_out)
%                     imagesc(im_edge(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                
                    
                subtightplot(2,4,6)
                    hold all
                    
                    imagesc(im_close(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(2,4,7)
                    hold all
                    
                    imagesc(im_fill(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
                subtightplot(2,4,8)
                    hold all
                    
                    imagesc(im_open(:,:,frame_idx))
                    
                    colormap gray
                    
                    axis image
                    set(gca,'Ydir','Reverse')

                    set(gca,'Color','white')
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                    
%                     error
                    pause
        end
        
        multiWaitbar('Segmenting...',frame_idx/num_frames);
    end
    
    im_mask = im_open;
    
    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% infer the surrounding bounding box for all wells put together
    
    multiWaitbar('CloseAll');
    multiWaitbar('Detecting objects...',0);
    
    objects = cell(1,num_frames);
    
    %%% object detection to find extrema of well spacing in this frame
    for frame_idx = 1:num_frames
        im_seg(:,:,frame_idx) = im_bf(:,:,frame_idx) .* im_mask(:,:,frame_idx);
        
        L = bwlabeln(im_seg(:,:,frame_idx));
        objects{frame_idx} = regionprops(L,im_seg(:,:,frame_idx),'Extrema','Area','Centroid');
        
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
    
    mask_final = zeros(size(im_bf));
    im_seg_final = zeros(size(im));
    
    %%% repair extrema matrices for frames that didn't segment well by
    %%% borrowing extrema from neighboring frames
    
    % sanity check
    
    extrema_extrema_x(extrema_extrema_x(:,1) + well_array_width > size(im_bf,2),:) = NaN;
    extrema_extrema_y(extrema_extrema_y(:,1) + well_array_height > size(im_bf,1),:) = NaN;
    
    mask_width = extrema_extrema_x(:,2) - extrema_extrema_x(:,1);
    mask_height = extrema_extrema_y(:,2) - extrema_extrema_y(:,1);
    
    extrema_extrema_x(mask_width < 0.9*well_array_width | mask_width > 1.1*well_array_width,:) = NaN;
    extrema_extrema_y(mask_height < 0.9*well_array_height | mask_height > 1.1*well_array_height,:) = NaN;
    
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
        
        for row_idx = 1:options.well_counts(1)
            for col_idx = 1:options.well_counts(2)
                start_i = round(extrema_extrema_y(frame_idx,1)) + (row_idx-1) * well_well_dist_height;
                end_i = start_i + options.well_height;
                start_j = round(extrema_extrema_x(frame_idx,1)) + (col_idx-1) * well_well_dist_width + options.processing_options.wseg_otsu2_left_offset;
                end_j = start_j + options.well_width;

                mask_final(start_i:end_i, start_j:end_j, frame_idx) = 1;
            end 
        end
        
        for channel_idx = 1:num_channels
            im_seg_final(:,:,frame_idx,channel_idx) = im(:,:,frame_idx,channel_idx) .* mask_final(:,:,frame_idx);
        end
    end
    
    % detect actual wells based on our inferred mask
    
    for frame_idx = 1:num_frames
        L = bwlabeln(mask_final(:,:,frame_idx));
        objects = regionprops(L,'Centroid','Extrema');
        
        if numel(objects) ~= prod(options.well_counts)
            error('The segmentation did not work properly...')
        end
        
        well_segmentation_results_struct.frame(frame_idx).all_objects = objects;
        well_segmentation_results_struct.frame(frame_idx).good_objects = objects;
    end
    
%     well_segmentation_results_struct.im_seg_final = im_seg_final;
end
