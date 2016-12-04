function [ output_args ] = test_radial_sym( input_args )
%TEST_RADIAL_SYM Summary of this function goes here
%   Detailed explanation goes here

% im_all = imread('single_molecule_signal.tif');
% im_all = imread('hard_brightfield.tif');
% im_all = imread('em_watershed_test.tif');
% im_all = imread('signal4.tif');

im_all = zloadim('em_watershed_test.tif');

im = im_all(:,:,1);
im = mat2gray(im);

processopt = 'spatialfilter';
processparam = [1 10];
thresh = 0.99;
fitstr = {'radial','none'};
try1pernhood = 0;
pad = 2*processparam(2);

im_pad = padarray(im,[pad pad]);

im_bpass = bpass(im_pad,processparam(1),processparam(2));

im_bpass = im_bpass(pad+1:end-pad,pad+1:end-pad);

im_fft = mat2gray(log(abs(fftshift(fft2(im)))+1));
im_bpass_fft = mat2gray(log(abs(fftshift(fft2(im_bpass)))+1));

% [cdf,x_loc] = ecdf(im_bpass(:));
    
frac = numel(im_bpass(im_bpass==0)) / numel(im_bpass(:));

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
        imagesc(im_bpass==0)
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

