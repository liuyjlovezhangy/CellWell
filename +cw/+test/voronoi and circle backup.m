function test_em_watershedding
%     im = imread('C3-Timelapse01.ome.tiff kept stack.tif');
%     im = imread('C2-Timelapse01.ome.tiff kept stack.tif');
%     im = imread('em_watershed_test.tif');
%     im = imread('joined_cells.tif');
    im = imread('joined_cells_multi.tif');
    
    im = mat2gray(im);

    % circle detection
    
%     h = fspecial('gaussian',[10 10], 100);
%     im_expanded = imfilter(im,h, 'replicate');
    im_expanded = im;

    im_expanded = imadjust(im_expanded,[],[],0.50);

    im_expanded = imresize(im_expanded,2);
    
%     h = fspecial('gaussian',[10 10], 100);
%     im_expanded = imfilter(im_expanded,h, 'replicate');
    
    im_expanded = convolveGaussian(im_expanded,1);

    im_expanded = imsharpen(im_expanded);

    [centers,radii] = imfindcircles(im_expanded,[5 20],'Sensitivity',0.8,'EdgeThreshold',0.1);
    
    % generate first pass mask
    
    cseg_tracking_params.threshold_density = 0.1;
    cseg_tracking_params.peak_stringency = 'low';
    cseg_tracking_params.threshold_smoothing = 'off';
    
    level = cw.process.adaptiveThresholdFinder(im_expanded, 1.5, cseg_tracking_params, []);
    im_thresh = threshold3D(im_expanded, level);
    
    % kill edge effect
    
    im_thresh(1,:) = 0;
    im_thresh(end,:) = 0;
    im_thresh(:,1) = 0;
    im_thresh(:,end) = 0;
    
    se = strel('disk',6);
    im_thresh = imopen(im_thresh,se);
%     im_thresh = imclearborder(im_thresh);

    % Clear very small objects
    im_thresh = bwareaopen(im_thresh, 100);
    
%     im_thresh = imresize(im_thresh,0.5);
    
    %%%
    
    % remove circles not in mask
    
%     centers = centers / 2;
    
    L = bwlabeln(im_thresh,4);
    objects = regionprops(L,'ConvexHull','Extrema','Image');
    
    circle_ids = NaN * ones(1,size(centers,1));
    
    contained = zeros(1,size(centers,1));
    
    im_thresh_final = zeros(size(im_thresh));
    
    for object_idx = 1:numel(objects)
        ch = objects(object_idx).ConvexHull;
        
        temp_contained = inpolygon(centers(:,1),centers(:,2),ch(:,1),ch(:,2))';
        
        circle_ids(temp_contained) = object_idx;
        
        contained = contained | temp_contained;
        
        temp_centers_contained = centers(temp_contained,:);
        
        %%%%%%%%%%%%
        % make an image the size of the original with just this object
        
        temp_image = zeros(size(im_thresh));
            
        start_i = min(ceil(objects(object_idx).Extrema(:,2)));
        start_j = min(ceil(objects(object_idx).Extrema(:,1)));

        imloc_i = [start_i:start_i + size(objects(object_idx).Image,1)-1];
        imloc_j = [start_j:start_j + size(objects(object_idx).Image,2)-1];

        temp_image(imloc_i,imloc_j) = objects(object_idx).Image;
        
        %%%%%%%%%%%%
        
        if sum(temp_contained) == 2

            % find where to put the line perpendicular to the two points
            
            line_center = mean(temp_centers_contained,1);
            
            m = diff(temp_centers_contained(:,2)) / diff(temp_centers_contained(:,1));
            orth_m = -1/m;
            
            linefun = @(x) orth_m * x - orth_m * line_center(1) + line_center(2);
            
            j_idcs=round(linspace(1,size(temp_image,2),sqrt(size(temp_image,1)^2 + size(temp_image,2)^2))); % the line will never need to be bigger than the image diagonal
            i_idcs = round(linefun(j_idcs));
            
            for idx = 1:numel(i_idcs)
                if i_idcs(idx) > size(temp_image,1)
                    continue
                end
                
                if i_idcs(idx) < 1
                    continue
                end
                
                temp_image(i_idcs(idx),j_idcs(idx)) = 0;
            end
            
%             se = strel('disk',4);
%             temp_image = imopen(temp_image,se);
        else
            
            temp_centers_contained = [temp_centers_contained(:,1) - start_j,temp_centers_contained(:,2) - start_i];
            
            ti = voronoi2mask(temp_centers_contained(:,1),temp_centers_contained(:,2),size(objects(object_idx).Image));

            [vxx,vy] = voronoi(temp_centers_contained(:,1),temp_centers_contained(:,2));
            
            masked_ti = ti;
            masked_ti(imcomplement(ti & objects(object_idx).Image)) = 0;

            
            figure(2314)
            clf
            
            subtightplot(1,3,1)
                hold all

                imagesc(objects(object_idx).Image)

                for center_idx = 1:size(temp_centers_contained,1)
                    plot(temp_centers_contained(center_idx,1),temp_centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                end

                axis image
                set(gca,'Ydir','Reverse')
                axis off 

                colormap gray
              
            subtightplot(1,3,2)
                hold all

                imagesc(ti)

                for center_idx = 1:size(temp_centers_contained,1)
                    plot(temp_centers_contained(center_idx,1),temp_centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                end

                axis image
                set(gca,'Ydir','Reverse')
                axis off 

                colormap gray
                
            subtightplot(1,3,3)
                hold all

                imagesc(masked_ti)

                for center_idx = 1:size(temp_centers_contained,1)
                    plot(temp_centers_contained(center_idx,1),temp_centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                end
                
                voronoi(temp_centers_contained(:,1),temp_centers_contained(:,2));

                axis image
                set(gca,'Ydir','Reverse')
                axis off 

                colormap gray
            
            return
        end
        
        im_thresh_final = im_thresh_final | temp_image;
    end
       
    centers_contained = centers(contained,:);
    
    L = bwlabeln(im_thresh_final,4);
    objects_final = regionprops(L,'ConvexHull','Extrema','Image');
    
    %%%

    figure(10394)
    clf
    
        subtightplot(2,3,1)
            hold all
            
            imagesc(im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,2)
            hold all
            
            imagesc(im_expanded)
            
            
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,3)
            hold all
            
            imagesc(im_expanded)
            viscircles(centers, radii,'EdgeColor','b');
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,4)
            hold all
            
            imagesc(im_thresh)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,5)
            hold all
            
            imagesc(im_thresh)
            
            for center_idx = 1:size(centers_contained,1)
                plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
            end
            
            for object_idx = 1:numel(objects)
                pts = objects(object_idx).ConvexHull;
                
                plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
            end
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,6)
            hold all
            
            imagesc(im_thresh_final)
            
            for center_idx = 1:size(centers_contained,1)
                plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
            end
            
            for object_idx = 1:numel(objects_final)
                pts = objects_final(object_idx).ConvexHull;
                
                plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
            end
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
end

function argh2    
    centers_process = round(fliplr(centers));
    
    centers_im = zeros(size(im_thresh));
    
    for center_idx = 1:size(centers,1)
        centers_im(sub2ind(size(centers_im),centers_process(center_idx,1),centers_process(center_idx,2))) = 1;
    end
    
    centers_im = double(centers_im);
%     h = fspecial('gaussian',[100 100], 100);
%     gauss_im = imfilter(centers_im,h, 'replicate');
    gauss_im = mat2gray(convolveGaussian(centers_im,5));
    %
    
%     cell_mask_em = imextendedmax(im_seg_water, 0.9);

    % D is the input to watershed after the initial
    % segmented image is cleaned up

    D_im = gauss_im;

    D_im = imcomplement(D_im);

%     D_im = imimposemin(D_im, ~cur_cell_mask | cell_mask_em);

    D_im(D_im > 0.9) = -inf;

    water_im = watershed(D_im);


    
    figure(10394)
    clf
    
        subtightplot(2,3,1)
            hold all
            
            imagesc(im_thresh)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,2)
            hold all
            
            imagesc(im_thresh)
            
            viscircles(centers, radii,'EdgeColor','b');
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,3)
            hold all
            
            imagesc(centers_im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,4)
            hold all
            
            imagesc(gauss_im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,5)
            hold all
            
            imagesc(D_im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,6)
            hold all
            
            L_im = water_im;
            L_im(water_im == 1) = 0;
            L_im(water_im == 0) = 100;
            
            L_im = label2rgb(L_im,'hsv',[0 0 0]);
            imagesc(L_im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap hsv
            
            freezeColors
end

function argh
    H_start = 0.1;
    H_final = 0.6;
    
    delta_H = 0.01;
    
    cseg_tracking_params.threshold_density = 0.1;
    cseg_tracking_params.peak_stringency = 'low';
    cseg_tracking_params.threshold_smoothing = 'off';
    
    %%%%% THRESHOLDING
    
    level = cw.process.adaptiveThresholdFinder(im, 1.5, cseg_tracking_params, []);
    im_thresh = threshold3D(im, level);
    
    %%%%% MORPHOLOGY
    
    % link together objects

    se = strel('disk',1);
    dt_mask = imclose(im_thresh,se);

    % Delete very small connections between cells

    se = strel('disk',1);
    dt_mask = imopen(dt_mask,se);

    % Clear very small objects
    dt_mask = bwareaopen(dt_mask, 10);

    %%%%% WATERSHEDDING
    
	D = bwdist(~dt_mask);
    D = mat2gray(D);
    
    H_vals = H_start:delta_H:H_final;
    
    N1s = zeros(1,numel(H_vals));
    N2s = N1s;
    
    for H_idx = 1:numel(H_vals)
    
        D_em = imextendedmax(D, H_vals(H_idx));

        L_cur = bwlabeln(D_em,4);
        N1s(H_idx) = max(L_cur(:));
    
        D_em_water = -D_em;
        D_em_water(~dt_mask) = -inf;
        L_water = watershed(D_em_water);
        
        N2s(H_idx) = max(L_water(:))-2;
    end
    
    [~,H_best_idx] = min(abs(N2s-N1s));
    H_best = H_vals(H_best_idx)
    
    D_em = D;imextendedmax(D, H_best);

    D_em_water = -D_em;
    D_em_water(~dt_mask) = -inf;
    L_water = watershed(D_em_water);
    
    figure(1342)
    clf
    
        hold all
        
        plot(H_vals, N1s)
        plot(H_vals, N2s)
        
        line([H_best, H_best],ylim)
        
        legend('N1','N2','Location','Best')
    
    
    %%%%% plot
    
    figure(223)
    clf
    
        subtightplot(2,3,1)
            hold all
            
            imagesc(im)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,2)
            hold all
            
            imagesc(im_thresh)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,3)
            hold all
            
            imagesc(dt_mask)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,4)
            hold all
            
            imagesc(D)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,5)
            hold all
            
            imagesc(D_em)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,6)
            hold all
            L_im = L_water;
            L_im(L_water == 1) = 0;
            L_im(L_water == 0) = 1;
            
            L_im = label2rgb(L_im,'hsv',[0 0 0]);
            imagesc(L_im)

            axis image
            set(gca,'Ydir','Reverse')
            axis off 

%             colormap hsv
    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
end