function im_rotate = rotate( im, angle_channel )
    if ~exist('angle_channel','var')
        angle_channel = 1;
    end

    noise_sz = 2;
    obj_sz = 50;
    pad = obj_sz * 3;
    
    angles = [];
    
    multiWaitbar('Detecting rotation angle...',0);
    
    num_frames = size(im,3);
    num_channels = size(im,4);
    
    for frame_idx = 1:min(num_frames,15)
        cur_im = mat2gray(wiener2(im(:,:,frame_idx,angle_channel)));
%         cur_im = mat2gray(im(:,:,frame_idx,angle_channel));
        
%         im_pad = padarray(cur_im,[pad pad],NaN);
%         im_bpass = bpass(im_pad,noise_sz,obj_sz);
%         cur_im = mat2gray(im_bpass(pad+1:end-pad,pad+1:end-pad));
        
%         figure(13423);clf;imagesc(im_bpass);colormap gray
%         pause
        
        cur_angle = horizon(cur_im,0.1,'hough');
        angles = [angles, cur_angle];
        
        multiWaitbar('Detecting rotation angle...',frame_idx/min(num_frames,15));
    end

    angle = mean(angles)
    
    multiWaitbar('CloseAll');
    multiWaitbar('Rotating...',0);
    
    im_rotate = zeros([size(imRotateCrop(im(:,:,1,angle_channel),-angle)), num_frames, num_channels]);
    
    for frame_idx = 1:num_frames
        for channel_idx = 1:num_channels
            im_rotate(:,:,frame_idx,channel_idx) = imRotateCrop(im(:,:,frame_idx,channel_idx),-angle);
        end
        
        multiWaitbar('Rotating...',frame_idx/num_frames);
    end
    multiWaitbar('CloseAll');
end

