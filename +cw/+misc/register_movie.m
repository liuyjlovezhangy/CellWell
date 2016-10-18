function register_movie(input_file,output_file,align_channel)
%     im = load_ome_tiff('movies/Timelapse01.ome.tiff');
    
    im = load_ome_tiff(input_file);

    if ~exist('align_channel','var')
        align_channel = 4;
    end
    
    usfac = 1;
    
    im_registered = zeros(size(im));
    
    im_ref = im(:,:,1,align_channel);
    
    for frame_idx = 1:size(im,3)
        frame_idx
        
        im_slice_cur = squeeze(im(:,:,frame_idx,:));
        
        [output, Greg] = dftregistration(fft2(im_ref),fft2(im_slice_cur(:,:,align_channel)),usfac);

        shiftI = round(output(3));
        shiftJ = round(output(4));
        
        im_registered(:,:,frame_idx,:) = shift_image(im_slice_cur,shiftI,shiftJ);
    end
    
%     figure(18342)
%     clf
% 
%         for frame_idx = 1:size(im,3)
%             subplot(1,2,1)
%                 cla
%                 imagesc(im(:,:,frame_idx,align_channel))
% 
%                 axis equal
%                 axis tight
%                 colormap gray
% 
%                 set(gca,'Ydir','Reverse')
%                 
%             subplot(1,2,2)
%                 cla
%                 imagesc(im_registered(:,:,frame_idx,align_channel))
% 
%                 axis equal
%                 axis tight
%                 colormap gray
% 
%                 set(gca,'Ydir','Reverse')
%             
%             drawnow
%         end

    im_registered = uint16(im_registered);
   
    delete(output_file)
    
%     bfsave(im_registered,output_file,'XYTCZ');

    t = Tiff(output_file,'w');
    
    tagstruct.ImageLength = size(im_registered,1);
    tagstruct.ImageWidth = size(im_registered,2);
    tagstruct.Photometric = Tiff.Photometric.RGB;
    tagstruct.BitsPerSample = 16;
    tagstruct.SamplesPerPixel = 4;
    tagstruct.RowsPerStrip = 16;
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    tagstruct.Software = 'MATLAB';
    tagstruct.ExtraSamples = Tiff.ExtraSamples.Unspecified;
    
    t.setTag('SubIFD',size(im,3));
    
    for frame_idx = 1:size(im,3)
        t.setTag(tagstruct);
        t.write(squeeze(im_registered(:,:,frame_idx,:)));
        t.writeDirectory();
    end
    
    t.close();
end

function im_shifted = shift_image(im,shiftI,shiftJ)
	% im is x,y,channel

    im_shifted = circshift(im,[shiftI,shiftJ,0]);    

    if shiftI > 0
        im_shifted(1:shiftI,:) = 0;
    elseif shiftI < 0
        im_shifted((size(im,1)+shiftI):size(im,1),:) = 0;
    end
    if shiftJ > 0
        im_shifted(:,1:shiftJ) = 0;
    elseif shiftJ < 0
        im_shifted(:,(size(im,2)+shiftJ):size(im,2)) = 0;
    end
end

