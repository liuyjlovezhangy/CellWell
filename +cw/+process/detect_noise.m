function [is_noise,im_filtered_list,im_filtered] = detect_noise( im, opts )
    
    % For this to work consistently, your images should be normalized across
    % the entire channel with something like mat2gray. The most important
    % part is that the maximum value is the maximum of the signal for the entire original image. You
    % should not normalize separate images / wells independently

    do_plot=0;
    gap = 0.1;

    if ~isfield(opts,'thresh_mean')
        opts.thresh_mean = 1.2;
    end
    if ~isfield(opts,'thresh_stdev')
        opts.thresh_stdev = 0.5;
    end
    if ~isfield(opts,'do_wiener')
        opts.do_wiener = 1;
    end
    if ~isfield(opts,'do_blur')
        opts.do_blur = 1;
    end
    if ~isfield(opts,'blur_sigma')
        opts.blur_sigma = 10;
    end
    if ~isfield(opts,'blur_hsize')
        opts.blur_hsize = [10 10];
    end
    if ~isfield(opts,'blur_hsize')
        opts.blur_hsize = [10 10];
    end
    if ~isfield(opts,'entropy_nsize')
        opts.entropy_nsize = ones(9);
    end
    
    plot_frame = 1;
    
    %%% This works on i,j,frame (3D) stacks
    
    if do_plot
        figure(1)
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
    
    if opts.do_wiener
        for frame_idx = 1:size(im,3)
            im(:,:,frame_idx) = wiener2(im(:,:,frame_idx),[5 5]);
        end
        
        edge_sz = 3;
        
        im = im(edge_sz:end-edge_sz,edge_sz:end-edge_sz,:);
    end
    
    if do_plot && opts.do_wiener
        subtightplot(2,3,3,gap)
            imagesc(im(:,:,plot_frame))
            colormap gray
            axis image
            set(gca,'Ydir','Reverse')
            axis off
            title('Wiener image')
    end
    
    if opts.do_blur
        h = fspecial('gaussian',opts.blur_hsize,opts.blur_sigma);
        im = imfilter(im,h);
        
        edge_sz = ceil(max(opts.blur_hsize) ./ 2);
        
        % edge effects
        im = im(edge_sz:end-edge_sz,edge_sz:end-edge_sz,:);
    end
    
    if do_plot && opts.do_blur
        subtightplot(2,3,4,gap)
            imagesc(im(:,:,plot_frame))
            colormap gray
            axis image
            set(gca,'Ydir','Reverse')
            axis off
            title('Gaussian blur image')
    end

    im_filtered = entropyfilt(im,opts.entropy_nsize);
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
    
    is_noise = allnans | ((mean_entropy < opts.thresh_mean) & (std_entropy < opts.thresh_stdev));
    is_noise = is_noise';
    
    if do_plot
        error
    end
end

