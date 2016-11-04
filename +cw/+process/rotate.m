function im_rotate = rotate( im, angle_channel )
    if ~exist('angle_channel','var')
        angle_channel = 1;
    end

    angles = [];
    
    multiWaitbar('Detecting rotation angle...',0);
    
    num_frames = size(im,3);
    num_channels = size(im,4);
    
    for frame_idx = 1:num_frames
        cur_angle = horizon(im(:,:,frame_idx,angle_channel),0.1,'hough');
        angles = [angles, cur_angle];
        
        multiWaitbar('Detecting rotation angle...',frame_idx/num_frames);
    end

    angle = mean(angles);
    
    multiWaitbar('CloseAll');
    multiWaitbar('Rotating...',0);
    
    im_rotate = zeros([size(imrotate(im(:,:,1,angle_channel),-angle)), num_frames, num_channels]);
    
    for frame_idx = 1:num_frames
        for channel_idx = 1:num_channels
            im_rotate(:,:,frame_idx,channel_idx) = imrotate(im(:,:,frame_idx,channel_idx),-angle);
        end
        
        multiWaitbar('Rotating...',frame_idx/num_frames);
    end
    
end

