function [ output_args ] = test_gap_closing( input_args )
%TEST_GAP_CLOSING Summary of this function goes here
%   Detailed explanation goes here

    molarray1 = [0 0 1 1; 0 0 2 1; 0 0 3 1; ...
        1 0 5 1; 1 0 6 1; 1 0 7 1; 1 0 8 1; ...
        2 0 10 1; 2 0 11 1; 2 0 12 1];
    
    molarray2 = [100 0 14 1; 100 0 15 2; 105 0 17 2; 105 0 18 2];
    
    molarray = [molarray1;molarray2];
    
    linkparams.searchrad = 10;
    linkparams.gap_close = 1;
    linkparams.mode = 'conservative';
    
    [trajectories,obj_ids] = simple_tracking( molarray, linkparams )
    
    for t_idx = 1:numel(trajectories)
        trajectories{t_idx}
    end
end

