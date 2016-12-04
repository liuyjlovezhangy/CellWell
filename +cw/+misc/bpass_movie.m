function [ output_args ] = test_radial_sym( input_args )
%TEST_RADIAL_SYM Summary of this function goes here
%   Detailed explanation goes here

% im_all = imread('single_molecule_signal.tif');
% im_all = imread('hard_brightfield.tif');
% im_all = imread('em_watershed_test.tif');
% im_all = imread('signal4.tif');

% im_all = zloadim('movies/01_CART1_target/Microwells_4color-03(2)_s2_corrected.tiff',1);
% im_all = zloadim('movies/04_CART1_only/Microwells_4color-03(18)_s18_corrected.tif',1);
im_all = zloadim('movies/Microwells_4color-04(10)_s10.ome.tiff',1);

processopt = 'spatialfilter';
processparam = [2 7];
fitstr = {'radial','none'};
try1pernhood = 0;
pad = 2*processparam(2);
thresh = 0.01;

for channel_idx = 1:4

    channel_idx
    
    new_im = zeros(size(im_all,1),size(im_all,2),size(im_all,3));
    
    for frame_idx = 1:size(im_all,3)
        im = im_all(:,:,frame_idx,channel_idx);
        im = mat2gray(im);

        im_pad = padarray(im,[pad pad]);

        im_bpass = bpass(im_pad,processparam(1),processparam(2));

        new_im(:,:,frame_idx) = im_bpass(pad+1:end-pad,pad+1:end-pad);
    end
    
    write_im(new_im,['im_out_' num2str(channel_idx) '.tif']);
end

write_im(mat2gray(im_all(:,:,:,5)),['im_out_' num2str(5) '.tif']);

    function write_im(image,filename)
        delete(filename)
        for cur_frame_idx = 1:size(image,3)
            imwrite(image(:,:,cur_frame_idx),filename,'WriteMode','append');
        end
    end

return

im_fft = mat2gray(log(abs(fftshift(fft2(im)))+1));
im_bpass_fft = mat2gray(log(abs(fftshift(fft2(im_bpass)))+1));

% [cdf,x_loc] = ecdf(im_bpass(:));

frac = numel(im_bpass(im_bpass>=thresh)) / numel(im_bpass(:));

figure(1342);clf
    subtightplot(1,3,1)
        hold all
        imagesc(im)
        axis image
        axis off
    subtightplot(1,3,2)
        hold all
        imagesc(im_bpass)
        axis image
        axis off
    subtightplot(1,3,3)
        hold all
        imagesc(im_bpass>=thresh)
        axis image
        axis off
        
    colormap gray
        
    suptitle(num2str(frac))
        
figure(2039);clf
    subtightplot(1,2,1)
        hold all
        imagesc(im_fft)
        axis image
        axis off
    subtightplot(1,2,2)
        hold all
        imagesc(im_bpass_fft)
        axis image
        axis off
%     subplot(2,4,4)
%         hold all
%         hist(im(:),1000)
%         xlim([0 1])
%     subplot(2,4,5)
%         hold all
%         hist(im_bpass(im_bpass~=0),1000)
%         xlim([0 1])
%     subplot(2,4,6)
%         hold all
%         plot(x_loc,cdf)
%         set(gca,'xscale','log')
        
    

colormap gray
end

