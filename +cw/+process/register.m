function im_register = register( im, registration_channel )
    if ~exist('registration_channel','var')
        registration_channel = 1;
    end

    usfac = 1;

    multiWaitbar('CloseAll');
    multiWaitbar('Registering...',0);

    im_register = zeros(size(im));
    
    num_frames = size(im,3);
    num_channels = size(im,4);
    
    for frame_idx = 1:num_frames
        [output, ~] = dftregistration(fft2(im(:,:,1,registration_channel)),fft2(im(:,:,frame_idx,registration_channel)),usfac);

        shiftI = round(output(3));
        shiftJ = round(output(4));

        for channel_idx = 1:num_channels
            im_register(:,:,frame_idx,channel_idx) = shift_image(im(:,:,frame_idx,channel_idx),shiftI,shiftJ);
        end

        multiWaitbar('Registering...',frame_idx/num_frames);
    end
    
    multiWaitbar('CloseAll');
end

