function new_im = noise_blaster( im_all, back_sz )

    N = 2*back_sz + 1;
    bx = zeros(1,N)  + 1/N;
    by = bx';

    new_im = zeros(size(im_all));

    for frame_idx = 1:size(im_all,3)
        
        im_slice = im_all(:,:,frame_idx);
        im_slice = mat2gray(im_slice);
        
        im_slice = wiener2(im_slice);

        if back_sz
        
            box1 = conv2(im_slice,bx,'same');
            box_final = conv2(box1,by,'same');
            
        else
            box_final = zeros(size(im_slice));
        end

        new_im(:,:,frame_idx) = im_slice - box_final;

    end
    
end