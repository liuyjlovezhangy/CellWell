function [signal_detection_results_struct, validation_images] = detect_noise( well_tracking_results_struct, options )
    
    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Detecting wells and frames that contain signal [cells]...',0);
    
    is_noise_matrix = zeros(num_wells,num_channels,num_frames); % well, channel, frame
    detection_images = cell(1,num_wells);

    im_well = well_tracking_results_struct.wells(1).im_well;
    
    validation_images = cell(1,num_channels);
    
    for channel_idx = 1:num_channels
        validation_images{channel_idx} = zeros(size(im_well,1)*options.well_counts(1),size(im_well,2)*options.well_counts(2),num_frames,3);
    end
    
    for well_idx = 1:numel(well_tracking_results_struct.wells)
        im_well = well_tracking_results_struct.wells(well_idx).im_well;

        % entropy calc chops off pixels on img border due to artifacting
        [~,~,sizing_im] = local( im_well(:,:,1,1), options.processing_options );
        cur_detection_image = zeros(size(sizing_im,1), size(sizing_im,2), num_frames, num_channels); 

        for channel_idx = 1:num_channels

            [temp_noise_matrix,~,channel_im_filtered] = local( im_well(:,:,:,channel_idx), options.processing_options );

            is_noise_matrix(well_idx,channel_idx,:) = temp_noise_matrix;
            cur_detection_image(:,:,:,channel_idx) = channel_im_filtered;
            
            % add this to the validation images cell so we can plot later
            % to see if user wants to keep this result
            
            if options.ask_me
                scale = 0.1;
                
                for frame_idx = 1:num_frames
                    if is_noise_matrix(well_idx,channel_idx,frame_idx)
                        valim = cat(3,im_well(:,:,frame_idx,channel_idx),zeros(size(im_well,1),size(im_well,2)),scale*ones(size(im_well,1),size(im_well,2)));
                    else
                        valim = cat(3,im_well(:,:,frame_idx,channel_idx),scale*ones(size(im_well,1),size(im_well,2)),zeros(size(im_well,1),size(im_well,2)));
                    end

                    [well_i,well_j] = ind2sub(options.well_counts,well_idx);

                    i = 1 + (well_i-1) * options.well_height;
                    i_end = well_i * options.well_height;

                    j = 1 + (well_j-1) * options.well_width;
                    j_end = well_j * options.well_width;

                    validation_images{channel_idx}(i:i_end,j:j_end,frame_idx,:) = valim;
                end
            end
        end

        detection_images{well_idx} = cur_detection_image;

        multiWaitbar('Detecting wells and frames that contain signal [cells]...',well_idx / numel(well_tracking_results_struct.wells));
    end

    signal_detection_results_struct.detection_images = detection_images;
    signal_detection_results_struct.is_noise_matrix = is_noise_matrix;
    signal_detection_results_struct.detection_opts = options;
    
    multiWaitbar('CloseAll');
    drawnow
end

function [is_noise,im_filtered_list,im_filtered] = local(im,opts)
    % For this to work consistently, your images should be normalized across
    % the entire channel with something like mat2gray. The most important
    % part is that the maximum value is the maximum of the signal for the entire original image. You
    % should not normalize separate images / wells independently

    do_plot=0;
    gap = 0.1;
    
    plot_frame = 1;
    
    %%% This works on i,j,frame (3D) stacks
    
    if do_plot
        figure(15341)
        clf

        subtightplot(2,3,1,gap)
            imagesc(im(:,:,plot_frame))
            colormap gray
            axis image
            set(gca,'Ydir','Reverse')
            axis off
            title('Original image')

        dist = im(:,:,plot_frame);
        dist = dist(:);

        subtightplot(2,3,2,gap)
            hist(dist,20)
            title('Intensity histogram')
    end
    
    if opts.ndec_do_wiener
        for frame_idx = 1:size(im,3)
            im(:,:,frame_idx) = wiener2(im(:,:,frame_idx),[5 5]);
        end
        
%         edge_sz = 3;
%         
%         im = im(edge_sz:end-edge_sz,edge_sz:end-edge_sz,:);
    end
    
    if do_plot && opts.ndec_do_wiener
        subtightplot(2,3,3,gap)
            imagesc(im(:,:,plot_frame))
            colormap gray
            axis image
            set(gca,'Ydir','Reverse')
            axis off
            title('Wiener image')
    end
    
    if opts.ndec_do_blur
        h = fspecial('gaussian',opts.ndec_blur_hsize,opts.ndec_blur_sigma);
        im = imfilter(im,h);
        
        edge_sz = ceil(max(opts.ndec_blur_hsize) ./ 2);
        
        % edge effects
        im = im(edge_sz:end-edge_sz,edge_sz:end-edge_sz,:);
    end
    
    if do_plot && opts.ndec_do_blur
        subtightplot(2,3,4,gap)
            imagesc(im(:,:,plot_frame))
            colormap gray
            axis image
            set(gca,'Ydir','Reverse')
            axis off
            title('Gaussian blur image')
    end

    im_filtered = entropyfilt(im,opts.ndec_entropy_nsize);
    im_filtered_list = reshape(im_filtered,[size(im_filtered,1) * size(im_filtered,2), size(im_filtered,3)])';
    
    if do_plot
            subtightplot(2,3,5,gap)
                imagesc(im_filtered(:,:,plot_frame))
                colormap gray
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                title('Entropy image')

            subtightplot(2,3,6,gap)
                dist = im_filtered(:,:,plot_frame);
                dist = dist(:);
                dist = dist(dist~=0);
                hist(dist,20)
                title('Entropy filter hist')

                mean_entropy = mean(dist)

                line([mean_entropy mean_entropy],ylim,'LineStyle','--','Color','r','LineWidth',3)


        set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');

        drawnow
    
    end
    
    im_filtered_list(im_filtered_list == 0) = NaN;
    
    allnans = all(isnan(im_filtered_list),2);
    mean_entropy = nanmean(im_filtered_list,2);
    std_entropy = nanstd(im_filtered_list,0,2);
    
    is_noise = allnans | ((mean_entropy < opts.ndec_thresh_mean) & (std_entropy < opts.ndec_thresh_stdev));
    is_noise = is_noise';
    
    if do_plot
        error
    end
end

