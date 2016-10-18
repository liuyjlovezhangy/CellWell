function batch_process( )
    base_dir = 'movies';

    file_listing = dir([base_dir '/*.ome.tiff']);
    
    for file_idx = 1:numel(file_listing)
        filename = file_listing(file_idx).name;
        
        full_filename = [base_dir '/' filename];
        
        try
            cw.process(full_filename)
        catch e
        end
    end
end

