function [trajectories,obj_ids] = simple_tracking( molarray, link_params )
%DTRACK_LINKFRAMES Takes the individual (x,y) points of the Gaussian peaks
%found using DaoStorm and reconstructs these points into trajectories using
%a simple search radius algorithm
% Input:
%   molarray structure: size==(# molecules, 3) with columns: x,y,frame 
%
%   link_params.searchrad = # of pixels for search radius
%
%   link_params.gap_close
%
% Output:
%
%
%   There are two operating modes:
%       link_params.mode = 'greedy': always connect a point to its nearest
%       neighbor (unless another point in the frame is closer to the
%       candidate) [NOT IMPLEMENTED YET]
%
%       link_params.mode = 'conservative': throw out any ambiguous linkages
%       where either one point could connect to two candidates in the next
%       frame or two points could connect to the same candidate

    num_frames = max(molarray(:,3));
    
    if link_params.gap_close
        gap_close_r = 2*link_params.searchrad;
    else
        gap_close_r = -inf;
    end
    
    molarray = [molarray zeros(size(molarray,1),1)] ;    
     
%     disp ('Constructing molarray...')
    track_num = 1 ;
    if strcmp(link_params.mode,'greedy')
        error('Greedy linking is not implemented yet')
    elseif strcmp(link_params.mode,'conservative')
%         h=waitbar(0,'Tracking objects...') ;
        % define x,y of all molecules in frame 1 as molecules in the
        % previous frame
        mol_ind = molarray(:,3) == 1 ;
        mol_pre = molarray(mol_ind,[1 2]);
        new_track_n = sum(mol_ind) ; % number of new tracks
        track_ind_pre = [1: new_track_n] ; % assign track indices
        molarray(mol_ind, 5) = track_ind_pre ;
        track_num = track_num + new_track_n ; % update track index counter              
         
        for frame_idx = 2:num_frames
            % define variable for x,y of all molecules in the currect frame 
            mol_ind = molarray(:,3) == frame_idx ;
            mol_cur = molarray(mol_ind,[1 2]);
            track_ind_cur = zeros(size(mol_cur,1),1) ; % track indices of molecules in the current frame
            for j = 1: size(mol_pre,1)
                dist_sq = sum((bsxfun(@minus,mol_cur,mol_pre(j,:))).^2, 2) ;
                candi_ind = dist_sq <= link_params.searchrad^2 ; % find candidate molecules
                
                if sum(candi_ind)==1 && track_ind_cur(candi_ind)==0 % If the molecule has only one candidate & the cadidate is not assigned yet
                    track_ind_cur(candi_ind)= track_ind_pre(j) ; % link them together
                    
                elseif sum(candi_ind)==1 && ~(track_ind_cur(candi_ind)==0) % If the molecule has only one candidate & the cadidate is already assigned to other track
                    track_ind_cur(candi_ind)= track_num ;  %% Create a new track index if more than one point links to the same candidate
                    track_num = track_num+1 ;
                    
                elseif sum(candi_ind) > 1
                    a = track_ind_cur(candi_ind);
                    i = 1;
                    found = 0;
                    
                    while i < size(dist_sq(candi_ind), 1)
                        if a(i) == 0 && found == 0
                            a(i) = track_ind_pre(j);
                            found = 1;
                        end
                        i = i +1;
                    end
                    
                else
                end            
            end
            % Create a new track index for all unassigned molecules
            new_track_n = sum(track_ind_cur==0) ;
            track_ind_cur(track_ind_cur==0)= [track_num:track_num+ new_track_n-1] ;
            track_num = track_num + new_track_n ;
            mol_pre = mol_cur ;
            track_ind_pre = track_ind_cur ;
            molarray(mol_ind, 5) = track_ind_cur ;
%             waitbar(frame/num_frames)
        end
%         close(h);
%         drawnow
        
        % compile trajectories

        trajectories = {};
        obj_ids = [];
        
        ignore_track_idcs = zeros(1,max(molarray(:,5)));
        
        track_idcs = unique(molarray(:,5));
        
        for track_idx = 1:numel(track_idcs)
            
            if ignore_track_idcs(track_idx)
                continue
            end

            [connected_track,cur_object_ids,start_frame,~,used_idcs] = find_track_connection(molarray,track_idcs(track_idx),gap_close_r);
            
            for ui_idx = 1:numel(used_idcs)
                ignore_track_idcs(used_idcs(ui_idx)) = 1;
            end
            
            tlen = size(connected_track,2);
            
            if ~isfield(link_params,'total_num_frames')
                total_num_frames = num_frames;
            else
                total_num_frames = link_params.total_num_frames;
            end
            
            if ~isfield(link_params,'min_track_len')
                min_track_len = 1;
            else
                min_track_len = link_params.min_track_len;
            end

            if size(connected_track,2) >= min_track_len
            
                % pad with NaNs on both sides

                connected_track = [NaN * ones(2,start_frame-1), connected_track, NaN * ones(2,total_num_frames - (start_frame+tlen-1))];
                cur_object_ids = [NaN * ones(1,start_frame-1), cur_object_ids, NaN * ones(1,total_num_frames - (start_frame+tlen-1))];

                trajectories = [trajectories, {connected_track}];
                obj_ids = [obj_ids; cur_object_ids];
            end

        end

    else
        error(['Unrecognized mode: ' link_params.mode]);
    end
    
end

function [connected_track,cur_object_ids,start_frame,end_frame,used_idcs] = find_track_connection(molarray,molecule_idx,gap_close_r)
    used_idcs = molecule_idx;

    molarray_subset = molarray(molarray(:,5) == molecule_idx,:);
    
    current_track = molarray_subset(:,1:2)';
    cur_object_ids = molarray_subset(:,4)';
    
    num_frames = max(molarray(:,3));
    
    start_frame = molarray_subset(1,3);
    end_frame = molarray_subset(end,3);
    
    connected_track = current_track;
    
    if isinf(gap_close_r)
        return
    end
    
    if end_frame ~= num_frames
        gap_frame = end_frame + 1;
        
        % find all subsequent tracks that start at gap_frame + 1
        
        for next_molecule_idx = molecule_idx+1:max((molarray(:,5)))
            molarray_subset_next = molarray(molarray(:,5) == next_molecule_idx,:);
            
            next_molecule_track = molarray_subset_next(:,1:2)';
            
            next_molecule_start_frame = molarray_subset_next(1,3);
            next_molecule_start_gap_frame = next_molecule_start_frame - 1;
            
            if gap_frame == next_molecule_start_gap_frame
                % see if these molecules are within the gap distance
                
                molecule_pos = current_track(:,end);
                next_molecule_pos = next_molecule_track(:,1);
                
                if eucl_dist(molecule_pos,next_molecule_pos) < gap_close_r
                    % these trajectories should be linked
                    
                    [recursive_connected_track,recursive_well_ids,~,end_frame,used_idcs_temp] = ...
                        find_track_connection(molarray,next_molecule_idx,gap_close_r);
                    
                    used_idcs = [used_idcs, used_idcs_temp];
                    
                    estimated_gap_pos = mean([molecule_pos,next_molecule_pos],2);
                    
                    connected_track = [connected_track,...
                        estimated_gap_pos,...
                        recursive_connected_track];
                    
                    cur_object_ids = [cur_object_ids, NaN, recursive_well_ids];

                    break
                end
            end
        end
    end
end