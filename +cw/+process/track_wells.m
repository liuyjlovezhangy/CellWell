function well_tracking_results_struct = track_wells( im, well_segmentation_results_struct, o )
        
    num_frames = size(im,3);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Initial tracking

    % build a list of all well locations throughout time

    localization_array = [];

    for frame_idx = 1:num_frames
        for well_idx = 1:numel(well_segmentation_results_struct.frame(frame_idx).good_objects)
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

    % Get average position trajectory

    well_tracks_mat = tracks_to_matrix(well_tracks);
    well_tracks_mat_shifted = zeros(size(well_tracks_mat));

    starting_average_position = squeeze(mean(well_tracks_mat,1));

    im_shifted = zeros(size(im));

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

    wells.track = [];

    % measure the average x,y length of the wells in the first frame

    well_half_widths = [];
    well_half_heights = [];

    for obj_idx = 1:numel(well_segmentation_results_struct.frame(1).good_objects)
        props = well_segmentation_results_struct.frame(1).good_objects(obj_idx);

        well_half_widths = [well_half_widths, ceil((max(props.Extrema(:,1)) - min(props.Extrema(:,1))) ./ 2)];
        well_half_heights = [well_half_heights, ceil((max(props.Extrema(:,2)) - min(props.Extrema(:,2))) ./ 2)];
    end

    well_half_widths = mean(well_half_widths);
    well_half_heights = mean(well_half_heights);

    new_mask = zeros(size(im,1),size(im,2));

    wells = repmat(wells,[1,total_wells]);

    for well_idx = 1:total_wells
        wells(well_idx).track = well_tracks{well_idx};
        wells(well_idx).track_shifted = well_tracks_shifted{well_idx};

        wells(well_idx).left_boundary = mean(well_tracks_shifted{well_idx}(1,:)) - well_half_widths - o.wseg_extra_border_x;
        wells(well_idx).right_boundary = mean(well_tracks_shifted{well_idx}(1,:)) + well_half_widths + o.wseg_extra_border_x;
        wells(well_idx).bottom_boundary = mean(well_tracks_shifted{well_idx}(2,:)) - well_half_heights - o.wseg_extra_border_y;
        wells(well_idx).top_boundary = mean(well_tracks_shifted{well_idx}(2,:)) + well_half_heights + o.wseg_extra_border_y;

        % make new mask now based on well tracking

        i_idcs = floor(wells(well_idx).bottom_boundary:ceil(wells(well_idx).top_boundary));
        j_idcs = floor(wells(well_idx).left_boundary):ceil(wells(well_idx).right_boundary);

        new_mask(i_idcs,j_idcs) = 1;

        wells(well_idx).im_well = im_shifted(i_idcs,j_idcs,:,:);

    end

    well_tracking_results_struct.wells = wells;
    well_tracking_results_struct.mask = new_mask;
end

