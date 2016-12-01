function [ output_args ] = test_radial_sym( input_args )
%TEST_RADIAL_SYM Summary of this function goes here
%   Detailed explanation goes here

im_all = imread('single_molecule_signal.tif');
im = im_all(:,:,1);
im = mat2gray(im);
im_small = imresize(im,0.25);
f = frst2d(im_small, 5, 5, 0.5, 'bright');

figure;
subplot(2,1,1)
imshow(im_small,[])
subplot(2,1,2)
imshow(f,[])
end

