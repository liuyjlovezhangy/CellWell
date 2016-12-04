function cell_segmentation_results_struct = segment_cells_radial_sym( well_tracking_results_struct, signal_detection_results_struct, options )
    
    processopt = 'spatialfilter';
    processparam = [1.5 38];
    thresh = 0.001;
    fitstr = {'radial','none'};
    try1pernhood = 0;
    padsize = 30;

    warning('segment_cells_circle currently assumes the BF channel is the last.')

    warning('off','images:imfindcircles:warnForLargeRadiusRange')
    warning('off','images:imfindcircles:warnForSmallRadius')
    
    propts = options.processing_options;

    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    cell_masks = cell(1,num_wells);
    
    detected_cell_props = cell(num_wells,num_frames,num_channels);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Performing cell segmentation...',0);
    
    for well_idx = 8:num_wells
        
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        cur_well_im_thresh = zeros(size(cur_well_im,1),size(cur_well_im,2),num_frames,num_channels);

        multiWaitbar('Current well...',0);

        for frame_idx = 1:num_frames

            for channel_idx = options.nuclei_channel
                
                if signal_detection_results_struct.is_noise_matrix(well_idx,channel_idx,frame_idx)
                    im_thresh_final = zeros(size(cur_well_im,1),size(cur_well_im,2));
                else
                    
                    im_slice = cur_well_im(:,:,frame_idx,channel_idx);
                    
%                     figure(12342);clf;imagesc(im_slice);colormap gray;axis image
%                     im_slice_mod = wiener2(im_slice);
                    im_slice_mod = imresize(im_slice,2,'method','bicubic');
                    im_slice_mod = padarray(im_slice_mod,[padsize padsize]);
                    im_slice_mod = bpass(im_slice_mod,processparam(1),processparam(2));
                    im_slice_mod = convolveGaussian(im_slice_mod,3);

                    im_thresh = imdilate(imopen(im_slice_mod > 0.01,strel('disk',6)),strel('disk',4));

%                     im_slice_mod(:)
                    
%                     figure(1256512);clf;plot(ecdf(im_slice_mod(:)));%colormap gray;axis image

                    im_slice_thresholded = im_slice_mod;
                    im_slice_thresholded(im_thresh(:) == 0) = 0;

                    im_slice_thresholded = reshape(im_slice_thresholded,size(im_slice_mod));

                    objs= fo5_rp(im_slice_thresholded, processopt, processparam, thresh, fitstr, try1pernhood);
                    
                    im_slice_mod = im_slice_mod(padsize+1:end-padsize,padsize+1:end-padsize);
                    im_thresh = im_thresh(padsize+1:end-padsize,padsize+1:end-padsize);
                    im_slice_thresholded = im_slice_thresholded(padsize+1:end-padsize,padsize+1:end-padsize);
                    
%                     error
                    
                    if ~isempty(objs)
                        centers = objs(1:2,:)' - padsize;
                    else
                        centers = [];
                    end
                    
                    % now return thresh to its original size
                    
                    im_slice_mod = imresize(im_slice_mod,0.5);
                    im_thresh = imresize(im_thresh,0.5);
                    im_slice_thresholded = imresize(im_slice_thresholded,0.5);
                    
                    %%% remove circles not in mask

                    % scale centers and radii because of resize
                    
                    centers = centers / 2;
                    
                    region_map = bwlabeln(im_thresh,4);
                    objects = regionprops(region_map,'ConvexHull','Extrema','FilledImage');

                    circle_ids = NaN * ones(1,size(centers,1));

                    contained = zeros(1,size(centers,1));
                    
                    im_thresh_final = zeros(size(im_thresh));
                    
                    centers_contained = [];
                    
                    if numel(centers) ~= 0 && ~isempty(objects)
                        
                        for object_idx = 1:numel(objects)
                            ch = objects(object_idx).ConvexHull;

                            temp_contained = inpolygon(centers(:,1),centers(:,2),ch(:,1),ch(:,2))';

                            circle_ids(temp_contained) = object_idx;

                            contained = contained | temp_contained;

                            temp_centers_contained = centers(temp_contained,:);

                            %%%%%%%%%%%%%%%
                            % find location we'll paste modified image into

                            temp_image = zeros(size(im_thresh));

                            start_i = min(ceil(objects(object_idx).Extrema(:,2)));
                            start_j = min(ceil(objects(object_idx).Extrema(:,1)));

                            imloc_i = [start_i:start_i + size(objects(object_idx).FilledImage,1)-1];
                            imloc_j = [start_j:start_j + size(objects(object_idx).FilledImage,2)-1];

                            %%%%%%%%%%%%
                            % Segment the current image if there are two or
                            % more circles present

                            if sum(temp_contained) < 2
                                temp_image(imloc_i,imloc_j) = objects(object_idx).FilledImage;

                            elseif sum(temp_contained) == 2

                                ti = objects(object_idx).FilledImage;

                                % find where to put the line perpendicular to the two points

                                temp_centers_contained_shifted = [temp_centers_contained(:,1) - start_j, temp_centers_contained(:,2) - start_i];

                                line_center = mean(temp_centers_contained_shifted,1);

                                m = diff(temp_centers_contained_shifted(:,2)) / diff(temp_centers_contained_shifted(:,1));
                                orth_m = -1/m;

                                ifun = @(j) orth_m * j - orth_m * line_center(1) + line_center(2);
                                jfun = @(i) (i - line_center(2))/orth_m + line_center(1);

                                % at really steep/shallow angles you have
                                % to change the independent variable or you
                                % won't have enough points to slice the
                                % objects apart
                                
                                if abs(orth_m) < 1 % shallow angle
                                    i_start = ifun(1);
                                    i_end = ifun(size(ti,1));

                                    i_idcs = linspace(i_start,i_end,2*sqrt(size(ti,1)^2 + size(ti,2)^2));
                                    j_idcs = round(jfun(i_idcs));
                                    i_idcs = round(i_idcs);
                                else
                                    j_start = jfun(1);
                                    j_end = jfun(size(ti,2));

                                    j_idcs = linspace(j_start,j_end,2*sqrt(size(ti,1)^2 + size(ti,2)^2));
                                    i_idcs = round(ifun(j_idcs));
                                    j_idcs = round(j_idcs);
                                end

                                for idx = 1:numel(i_idcs)
                                    if (i_idcs(idx) > size(ti,1)) || (j_idcs(idx) > size(ti,2))
                                        continue
                                    end

                                    if (i_idcs(idx) < 1) || (j_idcs(idx) < 1)
                                        continue
                                    end

                                    ti(i_idcs(idx),j_idcs(idx)) = 0;
                                end

                                if propts.cseg_debug
                                    figure(13542)
                                    clf
                                        hold all

                                        imagesc(ti)

                                        plot(temp_centers_contained_shifted(:,1),temp_centers_contained_shifted(:,2),'og-','MarkerSize',30,'LineWidth',5)
                                        plot(line_center(1),line_center(2),'xm','MarkerSize',30,'LineWidth',5)

                                        axis image
                                        set(gca,'Ydir','Reverse')
                                        axis off 

                                        colormap gray

%                                         pause
                                end

                                se = strel('disk',2);
                                ti = imopen(ti,se);

                                temp_image(imloc_i,imloc_j) = ti;

                            elseif sum(temp_contained) > 2

                                temp_centers_contained = [temp_centers_contained(:,1) - start_j,temp_centers_contained(:,2) - start_i];

                                ti = voronoi2mask(temp_centers_contained(:,1),temp_centers_contained(:,2),size(objects(object_idx).FilledImage));

                                masked_ti = ti;
                                masked_ti(imcomplement(ti & objects(object_idx).FilledImage)) = 0;

                                se = strel('disk',2);
                                masked_ti = imopen(masked_ti,se);

                                temp_image(imloc_i,imloc_j) = logical(masked_ti);
                            end

                            im_thresh_final = im_thresh_final | temp_image;
                        end

                        centers_contained = centers(contained,:);
%                     else
%                         centers_contained = [];
                    end
                    
                    im_thresh_final = bwlabeln(im_thresh_final,4);
                    objects_final = regionprops(im_thresh_final,'ConvexHull','Centroid','FilledImage','BoundingBox','Extrema','Area','EquivDiameter');
                    
                    % clean up objects that do not meet minimum
                    % requirements
                    
                    to = [];
                    
                    for object_idx = 1:numel(objects_final)
                        if objects_final(object_idx).Area < propts.cseg_min_area
                            start_i = min(ceil(objects_final(object_idx).Extrema(:,2)));
                            start_j = min(ceil(objects_final(object_idx).Extrema(:,1)));

                            imloc_i = [start_i:start_i + size(objects_final(object_idx).FilledImage,1)-1];
                            imloc_j = [start_j:start_j + size(objects_final(object_idx).FilledImage,2)-1];

                            im_thresh_final(imloc_i,imloc_j) = 0;
                        else
                            to = [to, objects_final(object_idx)];
                        end
                    end
                    
                    objects_final = to;
                    
                    detected_cell_props{well_idx,frame_idx,channel_idx} = objects_final;
                    
                    %%%%%%%%%%%%%%%
                    % add a border to all the objects for visualization
                    
                    im_thresh_temp = im_thresh_final;
                    im_thresh_temp(im_thresh_temp(:) ~= 0) = im_thresh_temp(im_thresh_temp(:) ~= 0) + 1;
                    
                    im_thresh_final = reshape(im_thresh_temp,size(im_thresh_final));
                    
                    for object_idx = 1:numel(objects_final)

                        start_i = min(ceil(objects_final(object_idx).Extrema(:,2)));
                        start_j = min(ceil(objects_final(object_idx).Extrema(:,1)));

                        imloc_i = [start_i:start_i + size(objects_final(object_idx).FilledImage,1)-1];
                        imloc_j = [start_j:start_j + size(objects_final(object_idx).FilledImage,2)-1];
                        
                        border_image = bwperim(objects_final(object_idx).FilledImage);
                        
                        ti = im_thresh_final(imloc_i,imloc_j);
                        ti2 = ti;
                        ti2(border_image == 1) = 1;
                        ti = reshape(ti2, size(ti));
                        
                        im_thresh_final(imloc_i,imloc_j) = ti;

                    end
                    
                    if propts.cseg_debug
                        figure(10394)
                        clf

                            subtightplot(2,3,1)
                                hold all

                                imagesc(im_slice)
                                
                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors

                            subtightplot(2,3,2)
                                hold all

                                imagesc(im_slice_mod)
                                
                                for center_idx = 1:size(centers_contained,1)
                                    plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'mo','MarkerSize',20,'LineWidth',2)
                                end

                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors

                            subtightplot(2,3,3)
                                hold all

                                imagesc(im_slice_thresholded)
                                
                                for center_idx = 1:size(centers_contained,1)
                                    plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'mo','MarkerSize',20,'LineWidth',2)
                                end

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
                                    plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'mo','MarkerSize',20,'LineWidth',2)
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

                    %             for center_idx = 1:size(centers_contained,1)
                    %                 plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'mo','MarkerSize',20,'LineWidth',2)
                    %             end

                                for object_idx = 1:numel(objects_final)
                                    pts = objects_final(object_idx).ConvexHull;

                                    plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
                                    plot(objects_final(object_idx).Centroid(1),objects_final(object_idx).Centroid(2),'mo','MarkerSize',20,'LineWidth',2)
                                end

                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors
                        
                        suptitle(['Well: ' num2str(well_idx) ' Frame: ' num2str(frame_idx) ' Channel: ' num2str(channel_idx)])
                                
                        pause
                    end
                end

                cur_well_im_thresh(:,:,frame_idx,channel_idx) = im_thresh_final;
                
            end

            multiWaitbar('Current well...',frame_idx / num_frames);
        end

        cell_masks{well_idx} = cur_well_im_thresh;

        multiWaitbar('Performing cell segmentation...',well_idx / num_wells);
    end

    cell_segmentation_results_struct.cell_masks = cell_masks;
    cell_segmentation_results_struct.detected_cell_props = detected_cell_props;
    cell_segmentation_results_struct.cell_segmentation_opts = propts;
    
    multiWaitbar('CloseAll');
    drawnow
end
