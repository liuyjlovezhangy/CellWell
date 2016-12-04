function [ output_args ] = test_radial_sym( input_args )
%TEST_RADIAL_SYM Summary of this function goes here
%   Detailed explanation goes here

% im_all = imread('single_molecule_signal.tif');
% im_all = imread('hard_brightfield.tif');
% im_all = imread('em_watershed_test.tif');
% im_all = imread('signal4.tif');

im_all = imread('HARD_CHANGING_BF.tif');

im = im_all(:,:,1);
im = mat2gray(im);

processopt = 'spatialfilter';
processparam = [1 31];
thresh = 0.99;
fitstr = {'radial','none'};
try1pernhood = 0;

res = bpass(im,processparam(1),processparam(2));

% objs= fo5_rp(im, processopt, processparam, thresh, fitstr, try1pernhood);

    

figure(1342);clf
subtightplot(2,1,1)
imagesc(im)
axis image
subtightplot(2,1,2)
imagesc(res)
axis image

colormap gray
end

