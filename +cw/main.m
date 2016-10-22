function main(options)
    
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
    mkdir([full_filename '__analysis_results/runs'])
    
    analysis_timestamp = strrep(datestr(now),':','-');    
    mkdir([full_filename '__analysis_results/runs/' analysis_timestamp])
    mkdir([full_filename '__analysis_results/runs/' analysis_timestamp '/movies'])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%
    %%% Run the processing on the movies
    %%%
    
    processing_options.filename = options.filename;
    
    if options.DO_PROCESSING
        processing_options
        cw.process(options);
    end
    
    %%% Copy all processing intermediates to runs/
    
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
    
    %%% Copy all analysis to runs/
    
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
    
    %%% Copy all movies to runs/
end