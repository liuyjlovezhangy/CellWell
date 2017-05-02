function well_tracking_results_struct = track_wells( im, well_segmentation_results_struct, options )

    num_frames = size(im,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initial tracking

    % build a list of all well locations throughout time

    localization_array = [];

    for frame_idx = 1:num_frames
        frame_idx
        for well_idx = 1:numel(well_segmentation_results_struct.frame(frame_idx).good_objects)
            well_idx
            objects = well_segmentation_results_struct.frame(frame_idx).good_objects(well_idx);

            localization_array = [localization_array; [objects.Centroid, frame_idx, well_idx]];

        end
    end
    
    % link the well positions frame-to-frame

    linkparams.mode = 'conservative';
    linkparams.searchrad = 20;
    linkparams.gap_close = 1;

    [well_tracks,well_ids] = cw.process.simple_tracking( localization_array, linkparams );
    
    well_tracking_results_struct.well_tracks = well_tracks;
    well_tracking_results_struct.well_ids = well_ids;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Correct for stage drift
    %%%%if stage not drifting, don't correct

    % Get average position trajectory

    well_tracks_mat = tracks_to_matrix(well_tracks);
    well_tracks_mat_shifted = zeros(size(well_tracks_mat));

    starting_average_position = squeeze(mean(well_tracks_mat,1));

    im_shifted = zeros(size(im));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%choose if you want shift wells because of drift
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    shift_choice = 1;  %%%don't want to shift
    
    if shift_choice == 1
    for frame_idx = 1:num_frames

        im_slice_cur = squeeze(im(:,:,frame_idx,:));

        % calculate shift

        shift = starting_average_position(frame_idx,:) - starting_average_position(1,:);
        im_shift = fliplr(shift);

        % shift image
        im_shifted(:,:,frame_idx,:) = shift_image(im_slice_cur,-im_shift(1),-im_shift(2));

        % correct tracks

        well_tracks_mat_shifted(:,frame_idx,:) = squeeze(well_tracks_mat(:,frame_idx,:)) - repmat(shift,[size(well_tracks_mat,1), 1]);
    end
    
    else
        for frame_idx = 1:num_frames
        im_shifted = im(:,:,frame_idx,:);
        well_tracks_mat_shifted(:,frame_idx,:) = well_tracks_mat(:,frame_idx,:);
        end
    end

    well_tracks_shifted = {};

    for track_idx = 1:size(well_tracks_mat_shifted,1)
        well_tracks_shifted = [well_tracks_shifted, {squeeze(well_tracks_mat_shifted(track_idx,:,:))'}];
    end

    well_tracking_results_struct.well_tracks_shifted = well_tracks_shifted;

    well_tracking_results_struct.well_tracks_mat = well_tracks_mat;
    well_tracking_results_struct.well_tracks_mat_shifted = well_tracks_mat_shifted;

    well_tracking_results_struct.im_shifted = im_shifted;

    % Add a wells entry to the structure where r_s.wells(well_idx) is
    % itself a structure containing the actual well's trajectory and
    % others

    total_wells = max(well_ids(:));

    if total_wells ~= prod(options.well_counts)
        total_wells
        prod_options.well_counts
        error('Well tracking: did not find enough wells during tracking')
    end
    
    wells.track = [];

    new_mask = zeros(size(im,1),size(im,2));

    wells = repmat(wells,[1,total_wells]);

    for well_idx = 1:total_wells
        wells(well_idx).track = well_tracks{well_idx};
        wells(well_idx).track_shifted = well_tracks_shifted{well_idx};

        %wells(well_idx).left_boundary = floor(mean(well_tracks_shifted{well_idx}(1,:)) - options.well_width/2);
        wells(well_idx).left_boundary = floor(mean(well_tracks_shifted{well_idx}(1,1)) - options.well_width/2);
        wells(well_idx).right_boundary = wells(well_idx).left_boundary + options.well_width - 1;
        %wells(well_idx).bottom_boundary = floor(mean(well_tracks_shifted{well_idx}(2,:)) - options.well_height/2);
        wells(well_idx).bottom_boundary = floor(mean(well_tracks_shifted{well_idx}(1,2)) - options.well_height/2);
        wells(well_idx).top_boundary = wells(well_idx).bottom_boundary + options.well_height - 1;

        % make well image stacks based on well tracking

        i_idcs = wells(well_idx).bottom_boundary:wells(well_idx).top_boundary;
        j_idcs = wells(well_idx).left_boundary:wells(well_idx).right_boundary;

%         del_idcs = (i_idcs < 1) | (j_idcs < 1);
        
        i_idcs(i_idcs < 1) = [];
        j_idcs(j_idcs < 1) = [];
        
        new_mask(i_idcs,j_idcs) = 1;

        wells(well_idx).im_well = im_shifted(i_idcs,j_idcs,:,:);
        
        figure(319872);clf
        subtightplot(1,2,1)
                    hold all
                    imagesc(im_shifted(:,:,1,end))
                    
%                     plot(wells(well_idx).track_shifted(1,1),wells(well_idx).track_shifted(2,1),'mo','MarkerSize',20,'LineWidth',4)
                    plot(j_idcs,i_idcs(1),'g-','MarkerSize',20,'LineWidth',4)
                    j_idcs
                    i_idcs
                    
%                     plot(mean(selected_centers(:,1)),mean(selected_centers(:,2)),'rx','MarkerSize',40,'LineWidth',4)
                    
                    colormap gray

                    axis tight
                    axis equal
                    axis off
                    set(gca,'Ydir','Reverse')  
        subtightplot(1,2,2)
                    hold all
                    imagesc(wells(well_idx).im_well(:,:,1,end))
%                     plot(selected_centers(:,1),selected_centers(:,2),'mo','MarkerSize',20,'LineWidth',4)
%                     plot(mean(selected_centers(:,1)),mean(selected_centers(:,2)),'rx','MarkerSize',40,'LineWidth',4)
                    
                    colormap gray

                    axis tight
                    axis equal
                    axis off
                    set(gca,'Ydir','Reverse')  
                    
                    error
    end

    well_tracking_results_struct.wells = wells;
    well_tracking_results_struct.mask = new_mask;
end

