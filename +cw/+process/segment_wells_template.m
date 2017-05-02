function [well_segmentation_results_struct, im_seg_final] = segment_wells_template(im, options)

    im_bf = im(:,:,:,options.bf_channel);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    well_well_dist_width = options.well_width + options.well_spacing_width;
    well_array_width = options.well_width * options.well_counts(2) + options.well_spacing_width * (options.well_counts(2)-1);
    
    well_well_dist_height = options.well_height + options.well_spacing_height;
    well_array_height = options.well_height * options.well_counts(1) + options.well_spacing_height * (options.well_counts(1)-1);
    
    template_filename = options.processing_options.wseg_template;
    
    template_im = zloadim(template_filename);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% init
    
    num_frames = size(im,3);
    num_channels = size(im,4);

    multiWaitbar('CloseAll');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% run segmentation

    mask_final = zeros(size(im_bf));
    im_seg_final = zeros(size(im));
    
    multiWaitbar('CloseAll');
    multiWaitbar('Segmenting...',0);
    
    center_avg_list = [];
    
    all_centers = {};
    
    for frame_idx = 1:num_frames
        
        cur_frame = im_bf(:,:,frame_idx);
        
        cur_frame = wiener2(cur_frame);

        I_NCCs = zeros([size(cur_frame), size(template_im,3)]);
        
        for template_idx = 1:size(template_im,3)
            template_slice = mat2gray(template_im(:,:,template_idx));
            
            %%% remove imageJ-added zero rows and columns
            
            rows = ~sum(template_slice,2);
            cols = ~sum(template_slice,1);
            
            template_slice(rows,:) = [];
            template_slice(:,cols) = [];
            
            template_slice = wiener2(template_slice);

            [~,I_NCC]=template_matching(template_slice,cur_frame);
            
            noise_sz = 1;
            obj_sz = 3;
            
            pad = 3*obj_sz;
             
            im_pad = padarray(I_NCC,[pad pad],NaN);
            im_bpass = bpass(im_pad,noise_sz,obj_sz);
            im_bpass = mat2gray(im_bpass(pad+1:end-pad,pad+1:end-pad));
            
            I_NCCs(:,:,template_idx) = im_bpass;
            
            if options.processing_options.wseg_debug
                
%                 figure(12342);clf;imagesc(template_slice);colormap gray;axis image;
%                 figure(12343);clf;imagesc(I_NCC);colormap gray;axis image;
%                 figure(12344);clf;imagesc(im_bpass);colormap gray;axis image;
%                 figure(12345);clf;imagesc(im_peaks);colormap gray;axis image;

%                 figure(13424);
%                     clf
%                     subtightplot(1,2,1)
%                         hold all
%                         imagesc(im_peaks)
%                         plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
%                         colormap gray
% 
%                         axis tight
%                         axis equal
%                         axis off
%                         set(gca,'Ydir','Reverse')
% 
%                     subtightplot(1,2,2)
%                         hold all
%                         imagesc(cur_frame)
%                         plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
%                         colormap gray
% 
%                         axis tight
%                         axis equal
%                         axis off
%                         set(gca,'Ydir','Reverse')   
            
%                 pause
            end
        end
        
        I_NCCs_sum = sum(I_NCCs,3);
        
        im_peaks = mat2gray(imclose(I_NCCs_sum,strel('disk',30)));
            
        centers=FastPeakFind(im_peaks,0);            
        centers = reshape(centers,2,[])';
        
        if size(centers,1) < options.well_counts(1)*options.well_counts(2)
            center_avg_list = [center_avg_list; [NaN NaN]];
        else
            idcs = sub2ind(size(im_peaks), centers(:,2), centers(:,1));
            values = im_peaks(idcs);
            [~,top_idcs] = sort(values,1,'descend');
            top_idcs = top_idcs(1:options.well_counts(1)*options.well_counts(2));
            selected_centers = centers(top_idcs,:);

            center_avg_list = [center_avg_list; mean(selected_centers,1)];
            
            all_centers = [all_centers, selected_centers];
%             disp test
        end
        
        if 0 && options.processing_options.wseg_debug
        
            figure(13424);
                clf
                subtightplot(1,2,1)
                    hold all
                    imagesc(im_peaks)
                    plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
                    plot(mean(selected_centers(:,1)),mean(selected_centers(:,2)),'rx','MarkerSize',40,'LineWidth',4)
                    
                    colormap gray

                    axis tight
                    axis equal
                    axis off
                    set(gca,'Ydir','Reverse')

                subtightplot(1,2,2)
                    hold all
                    imagesc(cur_frame)
                    plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
                    plot(mean(selected_centers(:,1)),mean(selected_centers(:,2)),'rx','MarkerSize',40,'LineWidth',4)
                    
                    colormap gray

                    axis tight
                    axis equal
                    axis off
                    set(gca,'Ydir','Reverse')  
            pause
        end
        
        multiWaitbar('Segmenting...',frame_idx/num_frames);
    end
    
    well_segmentation_results_struct = [];
    
    for frame_idx = 1:num_frames
        cur_centers = all_centers{frame_idx};
        
        objects = [];
        
        for obj_idx = 1:size(cur_centers,1)
            objects(obj_idx).Centroid = cur_centers(obj_idx,:);
        end
        
        well_segmentation_results_struct.frame(frame_idx).good_objects = objects;
    end
    
end
