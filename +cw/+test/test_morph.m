function test_morph()

    %%% Load image

%     im = im2double(imread('big_img.tif'));
    im = im2double(imread('C4-Timelapse02.ome.tif'));
    
    %%% Find otsu threshold levels and binarize
    
    num_levels = 2;
    
    [IDX,sep,idcs] = otsu(im,num_levels);

    im_mask_thresh = IDX == num_levels;
    
    %%% Image morphological operations to clean up
    
    im_mask_noborder = imclearborder(im_mask_thresh);
    
    se = strel('disk',10);
    im_mask_connect = imclose(im_mask_noborder,se);

    im_mask_noborder2 = imclearborder(im_mask_connect);
    
    im_mask_final = im_mask_noborder2;
    
    im_seg = im .* im_mask_final;
    
    %%% Object detection
    
    L = bwlabeln(im_seg);
    STATS = regionprops(L,im_seg,'Area','ConvexHull','MajorAxisLength','PixelList');
    
    %%% Plotting
    
    figure(81322)
    clf
    
    subplot(2,3,1)
        hold all
        
        imagesc(im)
        
        colormap gray
        axis tight
        axis equal
        set(gca,'Ydir','Reverse')
        
%         xlabel('x')
%         ylabel('y')
        title('Original image')
    
    subplot(2,3,2)
        hold all
        hist(im(:),2000)
        
        xlabel('Normalized pixel intensity')
        ylabel('# of pixels')
        title('Otsu thresholding')
        
        for idcs_idx = 1:numel(idcs)
            line([idcs(idcs_idx) idcs(idcs_idx)],ylim,'LineStyle','--','Color','r','LineWidth',3)
        end
        
        axis tight
        
    subplot(2,3,3)
        hold all
        
        imagesc(im_mask_thresh)
        
        colormap gray
        axis tight
        axis equal
        set(gca,'Ydir','Reverse')
        
%         xlabel('x')
%         ylabel('y')
        title('Otsu thresholded')
    
    subplot(2,3,4)
        hold all
        
        imagesc(im_mask_final)
        
        colormap gray
        axis tight
        axis equal
        set(gca,'Ydir','Reverse')
        
%         xlabel('x')
%         ylabel('y')
        title('Final mask')
        
    subplot(2,3,5)
        hold all
        
        imagesc(im)
        
        for obj_idx = 1:numel(STATS)

            pts = STATS(obj_idx).ConvexHull;

            plot(pts(:,1),pts(:,2),'g','LineWidth',2)

        end
        
        colormap gray
        axis tight
        axis equal
        set(gca,'Ydir','Reverse')
        
%         xlabel('x')
%         ylabel('y')
        title('Detected wells')
        
    subplot(2,3,6)
        hold all
        
        imagesc(im_seg)
        
        for obj_idx = 1:numel(STATS)

            pts = STATS(obj_idx).ConvexHull;

            plot(pts(:,1),pts(:,2),'g','LineWidth',2)

        end
        
        colormap gray
        axis tight
        axis equal
        set(gca,'Ydir','Reverse')
        
%         xlabel('x')
%         ylabel('y')
        title('Segmented wells')
        
    suptitle('Well detection performance')
        
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
%     set(gcf,'Position',[0 0 1400 700]);
%     fig = gcf;
%     fig.PaperPositionMode = 'auto';
end

