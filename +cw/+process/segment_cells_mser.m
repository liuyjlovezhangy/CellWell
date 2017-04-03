function [cell_segmentation_results_struct,validation_images] = segment_cells_point_source( well_tracking_results_struct, signal_detection_results_struct, options )
    
    warning('segment_cells_radial_sym currently assumes the BF channel is the last.')

    object_mask_overlap_threshold = 10;
    
    propts = options.processing_options;

    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    cell_masks = cell(1,num_wells);
    
    detected_cell_props = cell(num_wells,num_frames,num_channels);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Performing MSER cell segmentation...',0);
    
    validation_images = cell(1,num_channels * num_wells);
    
    for well_idx = 1:num_wells
        
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        cur_well_im_thresh = zeros(size(cur_well_im,1),size(cur_well_im,2),num_frames,num_channels);

        multiWaitbar('Current well...',0);

        for frame_idx = 1:num_frames

            for channel_idx = options.tracking_channels
                
                if signal_detection_results_struct.is_noise_matrix(well_idx,channel_idx,frame_idx)
                    im_thresh_final = zeros(size(cur_well_im,1),size(cur_well_im,2));
                else
                    
                    im_slice = cur_well_im(:,:,frame_idx,channel_idx);
                    
                    I = mat2gray(im_slice);
                    I = wiener2(I);
                    I = im2uint8(I);

                    [r,f] = vl_mser(I,'MinDiversity',0.4,'MaxVariation',0.2,'Delta',4,'BrightOnDark',1,'DarkOnBright',0) ;

                    % compute regions mask
                    M = zeros(size(I)) ;
                    for x=r'
                      s = vl_erfill(I,x) ;
                      M(s) = M(s) + 1;
                    end

                    M=imclose(M,strel('disk',3));
            %         M = imresize(M,1,'nearest');
%                     figure(013984);clf;imagesc(M);colormap gray

                    Mmask = M>0;
                    Mmask = bwareaopen(Mmask,14);
                    M(Mmask==0) = 0;

                    Mtemp = imcomplement(M);
                    Mtemp(M==0)=-inf;

                    L = watershed(Mtemp);
                    
                    cur_mask_water = L;
                    cur_corrected_mask = cur_mask_water;

                    cur_nowater_objects = regionprops(bwlabeln(Mmask,4),'Area','Centroid','PixelIdxList','ConvexHull','Image','BoundingBox','EquivDiameter','FilledImage');
                    cur_water_objects = regionprops(L,'Area','Centroid','PixelIdxList','ConvexHull','Image','BoundingBox','EquivDiameter','FilledImage');

                    % The largest area in this image will be the
                    % background unless something is very screwed up or
                    % you have loaded a carpet of cells

                    water_areas = [cur_water_objects.Area];
                    [~,background_idx] = max(water_areas);
                    cur_water_objects(background_idx) = [];

                    cur_objects_final = cur_water_objects;

                    for obj_nowater_idx = 1:numel(cur_nowater_objects)

                        % Can we find an object in the watershedded image that
                        % overlaps with the non-watershedded object?

                        % get pixel list of current nonwater object

                        nowater_pixels = cur_nowater_objects(obj_nowater_idx).PixelIdxList;

                        dont_add_nowater_obj = 0;

                        for obj_water_idx = 1:numel(cur_water_objects)
                            % put all checked watershedded objs into list

            %                 cur_objects_final = [cur_objects_final, cur_water_objects(obj_water_idx)];

                            % get pixel list of current water object

                            water_pixels = cur_water_objects(obj_water_idx).PixelIdxList;

                            % does the current nonwater object overlap with
                            % this water object by at least
                            % object_mask_overlap_threshold pixels?

                            object_overlap_area = sum(ismember(nowater_pixels,water_pixels));

                            if object_overlap_area > object_mask_overlap_threshold
                                % these are two apparently overlapping objects.
                                % this excluded nonwater from being added

                                dont_add_nowater_obj = 1;

                                break
                            end
                        end

                        if ~dont_add_nowater_obj
                            % This nonwater object did not correspond to a
                            % watershedded obj so we should added it to the
                            % final mask

                            next_obj_id = max(cur_corrected_mask(:))+1;

                            cur_corrected_mask(nowater_pixels) = next_obj_id;

                            % add a perimeter to the image

                            perim_im_size_matched = zeros(size(cur_corrected_mask));

                            perim_im = bwperim(cur_nowater_objects(obj_nowater_idx).Image);

                            bb = round(cur_nowater_objects(obj_nowater_idx).BoundingBox);

                            perim_im_size_matched(bb(2):bb(2)+bb(4)-1,bb(1):bb(1)+bb(3)-1) = perim_im;

                            cur_corrected_mask(logical(perim_im_size_matched)) = 0;
                            cur_objects_final = [cur_objects_final; cur_nowater_objects(obj_nowater_idx)];
                        end
                    end

                    final_map = cur_corrected_mask;

%                     cur_objects_final.Centroid
%                     error
                    
                    if ~isempty(cur_objects_final)
                        centers = [cur_objects_final.Centroid];
                    else
                        centers = [];
                    end
                    
                    radii = 5*ones(size(centers,1),1);
                    
                    if propts.cseg_debug
%                         figure(12302);clf;hold all;imagesc(imgLoG);colormap gray;axis image;set(gca,'Ydir','Reverse');axis off;set(gcf,'color','w')
%                         if ~isempty(centers)
%                             plot(centers(:,1),centers(:,2),'om','MarkerSize',20,'LineWidth',3)
%                         end
                    end
                    
                    %%% remove circles not in mask
                    
%                     region_map = bwlabeln(im_thresh,4);
%                     objects = regionprops(region_map,'ConvexHull','Extrema','FilledImage');
% 
%                     circle_ids = NaN * ones(1,size(centers,1));
% 
%                     contained = zeros(1,size(centers,1));
%                     
%                     im_thresh_final = zeros(size(im_thresh));
%                     
%                     centers_contained = [];
%                     
%                     if numel(centers) ~= 0 && ~isempty(objects)
%                         
%                         for object_idx = 1:numel(objects)
%                             ch = objects(object_idx).ConvexHull;
% 
%                             temp_contained = inpolygon(centers(:,1),centers(:,2),ch(:,1),ch(:,2))';
% 
%                             circle_ids(temp_contained) = object_idx;
% 
%                             contained = contained | temp_contained;
% 
%                             temp_centers_contained = centers(temp_contained,:);
% 
%                             %%%%%%%%%%%%%%%
%                             % find location we'll paste modified image into
% 
%                             temp_image = zeros(size(im_thresh));
% 
%                             start_i = min(ceil(objects(object_idx).Extrema(:,2)));
%                             start_j = min(ceil(objects(object_idx).Extrema(:,1)));
% 
%                             imloc_i = [start_i:start_i + size(objects(object_idx).FilledImage,1)-1];
%                             imloc_j = [start_j:start_j + size(objects(object_idx).FilledImage,2)-1];
% 
%                             %%%%%%%%%%%%
%                             % Segment the current image if there are two or
%                             % more circles present
% 
%                             if sum(temp_contained) < 2
%                                 temp_image(imloc_i,imloc_j) = objects(object_idx).FilledImage;
% 
%                             else
% 
%                                 ti = objects(object_idx).FilledImage;
% 
%                                 % find where to put the line perpendicular to the two points
% 
%                                 temp_centers_contained_shifted = [temp_centers_contained(:,1) - start_j, temp_centers_contained(:,2) - start_i];
% 
%                                 obj_map = zeros(size(ti));
%                                 
%                                 for center_idx = 1:size(temp_centers_contained_shifted,1)
%                                     obj_map(1+round(temp_centers_contained_shifted(center_idx,2)),1+round(temp_centers_contained_shifted(center_idx,1))) = 1;
%                                 end
% 
%                                 dist_map = bwdist(~obj_map);%%imcomplement(bwdist(obj_map));
%                                 
% %                                 if propts.cseg_debug
% %                                     figure(1622);clf;hold all;imagesc(dist_map);colormap gray;axis image;set(gca,'Ydir','Reverse');axis off;set(gcf,'color','w')
% %                                 end
%                                 
%                                 dist_map(~ti) = -inf;
%                                 dist_map = -dist_map;
%                                 
%                                 L_water = watershed(dist_map);
%                                 
%                                 if propts.cseg_debug
%                                     figure(1342);clf;imagesc(L_water);colormap gray;axis image;axis off;
%                                 end
%                                 
%                                 water_mask = L_water ~= 0;
%                                 
%                                 temp_image(imloc_i,imloc_j) = ti & water_mask;
% 
%                             end
% 
%                             im_thresh_final = im_thresh_final | temp_image;
%                         end
% 
%                         centers_contained = centers(contained,:);
% %                     else
% %                         centers_contained = [];
%                     end
%                     
%                     im_thresh_final = bwlabeln(im_thresh_final,4);
%                     objects_final = regionprops(im_thresh_final,'ConvexHull','Centroid','FilledImage','BoundingBox','Extrema','Area','EquivDiameter');
                    
                    % clean up objects that do not meet minimum
                    % requirements
                    
%                     to = [];
%                     
%                     for object_idx = 1:numel(objects_final)
%                         if objects_final(object_idx).Area < propts.cseg_min_area
%                             start_i = min(ceil(objects_final(object_idx).Extrema(:,2)));
%                             start_j = min(ceil(objects_final(object_idx).Extrema(:,1)));
% 
%                             imloc_i = [start_i:start_i + size(objects_final(object_idx).FilledImage,1)-1];
%                             imloc_j = [start_j:start_j + size(objects_final(object_idx).FilledImage,2)-1];
% 
%                             im_thresh_final(imloc_i,imloc_j) = 0;
%                         else
%                             to = [to, objects_final(object_idx)];
%                         end
%                     end
%                     
%                     objects_final = to;
                    
                    detected_cell_props{well_idx,frame_idx,channel_idx} = cur_objects_final;
                    
                    %%%%%%%%%%%%%%%
                    % add a border to all the objects for visualization
                    
%                     im_thresh_temp = im_thresh_final;
%                     im_thresh_temp(im_thresh_temp(:) ~= 0) = im_thresh_temp(im_thresh_temp(:) ~= 0) + 1;
%                     
%                     im_thresh_final = reshape(im_thresh_temp,size(im_thresh_final));
%                     
%                     for object_idx = 1:numel(objects_final)
% 
%                         start_i = min(ceil(objects_final(object_idx).Extrema(:,2)));
%                         start_j = min(ceil(objects_final(object_idx).Extrema(:,1)));
% 
%                         imloc_i = [start_i:start_i + size(objects_final(object_idx).FilledImage,1)-1];
%                         imloc_j = [start_j:start_j + size(objects_final(object_idx).FilledImage,2)-1];
%                         
%                         border_image = bwperim(objects_final(object_idx).FilledImage);
%                         
%                         ti = im_thresh_final(imloc_i,imloc_j);
%                         ti2 = ti;
%                         ti2(border_image == 1) = 1;
%                         ti = reshape(ti2, size(ti));
%                         
%                         im_thresh_final(imloc_i,imloc_j) = ti;

%                     end
%                     figure(1232);clf;imagesc(im_expanded);colormap gray;error;
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

%                             subtightplot(2,3,2)
%                                 hold all
% 
%                                 imagesc(im_expanded)
% 
%                                 axis image
%                                 set(gca,'Ydir','Reverse')
%                                 axis off 
% 
%                                 colormap gray

                                freezeColors

%                             subtightplot(2,3,3)
%                                 hold all
% 
%                                 imagesc(im_expanded_thresholded)
%                                 viscircles(centers, radii,'EdgeColor','b');
% 
%                                 axis image
%                                 set(gca,'Ydir','Reverse')
%                                 axis off 
% 
%                                 colormap gray
% 
%                                 freezeColors

                            subtightplot(2,3,4)
                                hold all

                                imagesc(final_map)

                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors

                            subtightplot(2,3,5)
                                hold all

                                imagesc(final_map)

                                for center_idx = 1:size(centers,1)
                                    plot(centers(center_idx,1),centers(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                                end

                                for object_idx = 1:numel(cur_objects_final)
                                    pts = cur_objects_final(object_idx).ConvexHull;

                                    plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
                                end

                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors

                            subtightplot(2,3,6)
                                hold all

                                imagesc(final_map)

                    %             for center_idx = 1:size(centers_contained,1)
                    %                 plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                    %             end

                                for object_idx = 1:numel(cur_objects_final)
                                    pts = cur_objects_final(object_idx).ConvexHull;

                                    plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
                                    plot(cur_objects_final(object_idx).Centroid(1),cur_objects_final(object_idx).Centroid(2),'xr','MarkerSize',10,'LineWidth',2)
                                end

                                axis image
                                set(gca,'Ydir','Reverse')
                                axis off 

                                colormap gray

                                freezeColors
                        
                        suptitle(['Well: ' num2str(well_idx) ' Frame: ' num2str(frame_idx) ' Channel: ' num2str(channel_idx)])
                                
                        set(gcf,'color','w')
                        
                        pause
                    end
                end

                cur_well_im_thresh(:,:,frame_idx,channel_idx) = final_map;
                
            end

            multiWaitbar('Current well...',frame_idx / num_frames);
        end

        cell_masks{well_idx} = cur_well_im_thresh;
        
%         all_threshold_levels{well_idx} = cur_threshold_levels;
%         all_thresh_xvals{well_idx} = cur_thresh_xvals;
%         all_thresh_yvals{well_idx} = cur_thresh_yvals;

        multiWaitbar('Performing cell segmentation...',well_idx / num_wells);
    end
    
    cell_segmentation_results_struct.cell_masks = cell_masks;
%     cell_segmentation_results_struct.threshold_levels = all_threshold_levels;
%     cell_segmentation_results_struct.thresh_xvals = all_thresh_xvals;
%     cell_segmentation_results_struct.thresh_yvals = all_thresh_yvals;
    cell_segmentation_results_struct.detected_cell_props = detected_cell_props;
    cell_segmentation_results_struct.cell_segmentation_opts = propts;
    
    multiWaitbar('CloseAll');
    drawnow

end
