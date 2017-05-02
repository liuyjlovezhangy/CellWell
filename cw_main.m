function cw_main(options_file)

%     if ~exist('options_file','var')
% 
%         [filename, path] = uigetfile('*.m', 'Select CellWell configuration script...');
%         
%         options_file = [path filename];
%         
%         if filename == 0
%             return
%         end
%     end
% 
%     try
%         asdf
%         options = evalf(options_file);
%     catch e
%        return 
%     end
    
    if ~exist('options_file','var')

        [filename, path] = uigetfile('*.mat', 'Select CellWell configuration .mat...');
        
        options_file = [path filename];
        
        if filename == 0
            return
        end
    end

    options = importdata(options_file);
    
    processing_options = options.processing_options;
    analysis_options = options.analysis_options;
    plot_options = options.plot_options;

    disp('CELL WELL ANALYSIS: RUNNING')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Directory setup
    %%%
    
    full_filename = options.filename;
    
    if ~exist(full_filename,'file')
        msgbox('ERROR: Could not find supplied movie file. Check current directory.')
        error('ERROR: Could not find supplied movie file. Check current directory.')
    end
    
    mkdir([full_filename '__analysis_results'])
    
    if options.NUKE_IT
        if options.NUKE_WARN
            response = questdlg('You have selected to NUKE the current output on this movie from all runs from orbit. Nuke it?');
        else
            response = 'Yes';
        end
        
        if strcmp(response,'Yes')
            disp('Nuking old output data from orbit...')
            
            rmdir([full_filename '__analysis_results'],'s')
        elseif strcmp(response,'No')
        else
            msgbox('ERROR: Set the NUKE_IT flag to 0 in the options.')
            error('ERROR: Set the NUKE_IT flag to 0 in the options.')
        end
    end
    
    mkdir([full_filename '__analysis_results/movies'])
%     mkdir([full_filename '__analysis_results/runs'])
%     
%     analysis_timestamp = strrep(datestr(now),':','-');    
%     mkdir([full_filename '__analysis_results/runs/' analysis_timestamp])
%     mkdir([full_filename '__analysis_results/runs/' analysis_timestamp '/movies'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Run the processing on the movies
    %%%
    
    if options.DO_PROCESSING
        
        processing_options.filename = options.filename;

        % convert start to the appropriate index from string

        sa = processing_options.start_at;

        switch sa
            case 'well segmentation'
                new_sa = 1;
            case 'well tracking'
                new_sa = 2;
            case 'noise detection'
                new_sa = 3;
            case 'cell segmentation'
                new_sa = 4;
            case 'cell tracking'
                new_sa = 5;
            case 'next'
                new_sa = 1e100; % inf:finite is not working in mats_to_delete indexing below...
            otherwise
                new_sa = 0;
        end
        options.processing_options.start_at = new_sa;

        % kill all intermediate .mats that are generated beyond the current
        % step we want to start at

        mats = {'input_movie.mat','well_segmentation.mat',...
            'well_tracking.mat','noise_detection.mat','cell_segmentation.mat',...
            'cell_tracking.mat','cell_interactions.mat'};

        mats_to_delete = mats((new_sa+1):5);

        for kill_idx = 1:numel(mats_to_delete)
            kill_file = [full_filename '__analysis_results/' mats_to_delete{kill_idx}];
            if exist(kill_file,'file')
                delete(kill_file)
            end
        end

        % do it
    
        processing_options
        cw.process(options);
    end
    
    %%% Copy all processing intermediates to runs/ folder
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Analyze cell data, calculate statistics...
    %%%
    
    if options.DO_ANALYSIS
        analysis_options
        cw.analyze(options);
    end
    
    %%% Copy all analysis to runs/ folder
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Plot figures
    %%%
    
    if options.DO_PLOTTING
        plot_options
        cw.plotting(options);
    end
    
    %%% Copy all movies to runs/ folder
end