function [cell_segmentation_results_struct,validation_images] = segment_cells_point_source( well_tracking_results_struct, signal_detection_results_struct, options )
    
    warning('segment_cells_radial_sym currently assumes the BF channel is the last.')

    noise_sz = 2;
    obj_sz = 10;
    
    propts = options.processing_options;

    num_wells = numel(well_tracking_results_struct.wells);
    num_frames = size(well_tracking_results_struct.wells(1).im_well,3);
    num_channels = size(well_tracking_results_struct.wells(1).im_well,4);
    
    cell_masks = cell(1,num_wells);
    
    detected_cell_props = cell(num_wells,num_frames,num_channels);
    
    all_threshold_levels = cell(1,num_wells);
    all_thresh_xvals = cell(1,num_wells);
    all_thresh_yvals = cell(1,num_wells);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Performing cell segmentation...',0);
    
    validation_images = cell(1,num_channels * num_wells);
    
    for well_idx = 1:num_wells
        
        cur_well_im = mat2gray(well_tracking_results_struct.wells(well_idx).im_well);
        cur_well_im_thresh = zeros(size(cur_well_im,1),size(cur_well_im,2),num_frames,num_channels);

        cur_threshold_levels = cell(num_frames,num_channels);

        cur_thresh_xvals = cell(num_frames,num_channels);
        cur_thresh_yvals = cell(num_frames,num_channels);

        multiWaitbar('Current well...',0);

        for frame_idx = 1:num_frames

            for channel_idx = options.tracking_channels
                
                if signal_detection_results_struct.is_noise_matrix(well_idx,channel_idx,frame_idx)
                    im_thresh_final = zeros(size(cur_well_im,1),size(cur_well_im,2));
                else
                    
                    im_slice = cur_well_im(:,:,frame_idx,channel_idx);
                    
                    noise_sz = 1;
                    obj_sz = 10;
    
%                     im_slice = imclose(im_slice,strel('disk',4));
                    
%                     im_slice = bpass(im_slice,noise_sz,obj_sz);
%                     im_slice = mat2gray(im_slice);
%                     im_slice(im_slice<0.01) = 0;
%                     im_slice = convolveGaussian(im_slice,3);

                    % set the threshold to be the kink in the histogram
                    % which separates bnkd from signal
                    % https://www.mathworks.com/matlabcentral/answers/57194-how-to-find-the-sharp-turn-in-a-2d-line-curve

%                     data = im_slice(:);
%                     num_bins = numel(data)^(1/2);
%                     bin_width = (max(data(:)) - min(data(:))) / num_bins;
%                     [hist_data,hist_loc] = histo(data,bin_width);
% 
%                     CDF = cumsum(hist_data); CDF = CDF / max(CDF);
%                     
%                     thresh = hist_loc(find(CDF>0.99,1));
                    
%                     figure(13083);clf;
%                         subplot(3,1,1);hold all;plot(hist_loc,hist_data)
%                         subplot(3,1,2);hold all;plot(hist_loc,CDF)
% %                         subplot(3,1,2);hold all;plot(hist_loc(1:end-1),mat2gray(hist_data_diff))
% %                         subplot(3,1,3);hold all;plot(hist_loc(1:end-2),mat2gray(hist_data_diff2))
%                         
%                     
%                     error
                    
%                     im_slice(im_slice<thresh) = 0;
% 
%                     im_thresh = im_slice > thresh;
%                     error
                    
%                     figure(2033);clf;hold all;imagesc(im_thresh);colormap gray;axis image;

%                     error
                    
                    sigma = 2;
                    
                    [pstruct, im_thresh, imgLM, imgLoG] = pointSourceDetection(im_slice, sigma, 'Alpha', 0.1,'WindowSize',2*sigma,...
                        'RedundancyRadius',5);%,'mask',im_thresh);
                    
                    imgLoG = mat2gray(imgLoG);
                    
                    imgLoG = imcomplement(imgLoG);
                    
                    imgLoG_ac = activecontour(imgLoG,imdilate(im_thresh,strel('disk',3)),20); 
%                     imgLoG_ac = imcomplement(imgLoG_ac);
                    
                    figure(03894);
                    clf
                    subplot(1,2,1);hold all;imagesc(imgLoG);colormap gray;axis image
                    subplot(1,2,2);hold all;imagesc(imgLoG_ac);colormap gray;axis image
                    
                    
%                     error
                    
                    im_slice(imgLM>0) = 0;    
                    imgLoG(imgLM>0) = 0;   
                    imgLM(imgLM>0) = 1;

                    comp_im = cat(3,imgLM,im_slice,zeros(size(im_slice)));
                    comp_im_log = cat(3,imgLM,imgLoG,zeros(size(im_slice)));

                    frame_idx
                    channel_idx
                    
                    figure(1033);clf;hold all;imagesc(im_slice);colormap gray;axis image;if ~isempty(pstruct);plot(pstruct.x,pstruct.y,'mo','MarkerSize',20);end;
                    figure(1034);clf;hold all;imagesc(im_thresh);colormap gray;axis image;if ~isempty(pstruct);plot(pstruct.x,pstruct.y,'mo','MarkerSize',20);end;
                    figure(1035);clf;hold all;imagesc(comp_im);colormap gray;axis image;if ~isempty(pstruct);plot(pstruct.x,pstruct.y,'mo','MarkerSize',20);end;
                    figure(1036);clf;hold all;imagesc(comp_im_log);colormap gray;axis image;if ~isempty(pstruct);plot(pstruct.x,pstruct.y,'mo','MarkerSize',20);end;
%                     
%                     pstruct
                    if ~isempty(pstruct)
                        centers = [pstruct.x', pstruct.y'];
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

                            else

                                ti = objects(object_idx).FilledImage;

                                % find where to put the line perpendicular to the two points

                                temp_centers_contained_shifted = [temp_centers_contained(:,1) - start_j, temp_centers_contained(:,2) - start_i];

                                obj_map = zeros(size(ti));
                                
                                for center_idx = 1:size(temp_centers_contained_shifted,1)
                                    obj_map(1+round(temp_centers_contained_shifted(center_idx,2)),1+round(temp_centers_contained_shifted(center_idx,1))) = 1;
                                end

                                dist_map = bwdist(~obj_map);%%imcomplement(bwdist(obj_map));
                                
%                                 if propts.cseg_debug
%                                     figure(1622);clf;hold all;imagesc(dist_map);colormap gray;axis image;set(gca,'Ydir','Reverse');axis off;set(gcf,'color','w')
%                                 end
                                
                                dist_map(~ti) = -inf;
                                dist_map = -dist_map;
                                
                                L_water = watershed(dist_map);
                                
                                if propts.cseg_debug
                                    figure(1342);clf;imagesc(L_water);colormap gray;axis image;axis off;
                                end
                                
                                water_mask = L_water ~= 0;
                                
                                temp_image(imloc_i,imloc_j) = ti & water_mask;

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

                    %             for center_idx = 1:size(centers_contained,1)
                    %                 plot(centers_contained(center_idx,1),centers_contained(center_idx,2),'xr','MarkerSize',10,'LineWidth',2)
                    %             end

                                for object_idx = 1:numel(objects_final)
                                    pts = objects_final(object_idx).ConvexHull;

                                    plot(pts(:,1),pts(:,2),'g-','LineWidth',3)
                                    plot(objects_final(object_idx).Centroid(1),objects_final(object_idx).Centroid(2),'xr','MarkerSize',10,'LineWidth',2)
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

                cur_well_im_thresh(:,:,frame_idx,channel_idx) = im_thresh_final;
                
            end

            multiWaitbar('Current well...',frame_idx / num_frames);
        end

        cell_masks{well_idx} = cur_well_im_thresh;
        
        all_threshold_levels{well_idx} = cur_threshold_levels;
        all_thresh_xvals{well_idx} = cur_thresh_xvals;
        all_thresh_yvals{well_idx} = cur_thresh_yvals;

        multiWaitbar('Performing cell segmentation...',well_idx / num_wells);
    end
    
    cell_segmentation_results_struct.cell_masks = cell_masks;
    cell_segmentation_results_struct.threshold_levels = all_threshold_levels;
    cell_segmentation_results_struct.thresh_xvals = all_thresh_xvals;
    cell_segmentation_results_struct.thresh_yvals = all_thresh_yvals;
    cell_segmentation_results_struct.detected_cell_props = detected_cell_props;
    cell_segmentation_results_struct.cell_segmentation_opts = propts;
    
    multiWaitbar('CloseAll');
    drawnow

end
