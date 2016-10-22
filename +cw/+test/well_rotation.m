function well_rotation
    I = imread('rotated_wells.tif');
    
    angle_fft = horizon(I,0.1);
    angle_h = horizon(I,0.1,'hough');
%     angle = horizon(I,1,'blot',100);
    
    I_rot_fft = imrotate(I,-angle_fft);
    I_rot_h = imrotate(I,-angle_h);
    
    num_otsu_levels = 2;
    
    [IDX,sep,well_segmentation_results_struct.otsu_idcs] = otsu(I,num_otsu_levels);
    im_mask_thresh = IDX == num_otsu_levels;
    
    [IDX,sep,well_segmentation_results_struct.otsu_idcs] = otsu(I_rot_fft,num_otsu_levels);
    im_mask_thresh_fft_rot = IDX == num_otsu_levels;
    
    [IDX,sep,well_segmentation_results_struct.otsu_idcs] = otsu(I_rot_h,num_otsu_levels);
    im_mask_thresh_h_rot = IDX == num_otsu_levels;
    
    figure(13042)
    clf
    
        subtightplot(1,3,1)
            hold all

            imagesc(I)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('pre rot')
            
        subtightplot(1,3,2)
            hold all

            imagesc(I_rot_fft)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title(['FFT rot ' num2str(angle_fft) ' degrees'])
            
        subtightplot(1,3,3)
            hold all

            imagesc(I_rot_h)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title(['Hough rot ' num2str(angle_h) ' degrees'])
                    
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
    
    figure(1512)
    clf
        subtightplot(2,3,1)
            hold all

            imagesc(I)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title('pre rot')
            
        subtightplot(2,3,2)
            hold all

            imagesc(I_rot_fft)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title(['FFT rot ' num2str(angle_fft) ' degrees'])
            
        subtightplot(2,3,3)
            hold all

            imagesc(I_rot_h)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
            title(['Hough rot ' num2str(angle_h) ' degrees'])
    
        subtightplot(2,3,4)
            hold all

            imagesc(im_mask_thresh)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
        
        subtightplot(2,3,5)
            hold all

            imagesc(im_mask_thresh_fft_rot)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
        subtightplot(2,3,6)
            hold all

            imagesc(im_mask_thresh_h_rot)
            colormap gray
            
            axis tight
            axis equal
            axis off
            set(gca,'Ydir','Reverse')
            
    set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
    set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
    set(gcf, 'color', 'white');
end

function based_on_object_properties
    clc;    % Clear the command window.
    close all;  % Close all figures (except those of imtool.)
    imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
    workspace;  % Make sure the workspace panel is showing.
    format long g;
    format compact;
    fontSize = 22;
    % Read in a demo image.
    folder = 'D:\Temporary Stuff';
    baseFileName = 'can_end.jpg';
    % Get the full filename, with path prepended.
    fullFileName = fullfile(folder, baseFileName);
    if ~exist(fullFileName, 'file')
        % Didn't find it there.  Check the search path for it.
        fullFileName = baseFileName; % No path this time.
        if ~exist(fullFileName, 'file')
            % Still didn't find it.  Alert user.
            errorMessage = sprintf('Error: %s does not exist.', fullFileName);
            uiwait(warndlg(errorMessage));
            return;
        end
    end
    rgbImage = imread(fullFileName);
    % Get the dimensions of the image.  numberOfColorBands should be = 3.
    [rows, columns, numberOfColorBands] = size(rgbImage);
    if numberOfColorBands > 1
        grayImage = rgb2gray(rgbImage);
    else
        grayImage = rgbImage;
    end
            % Display the original color image.
            subplot(3, 3, 1);
            imshow(grayImage);
            title('Original Image', 'FontSize', fontSize);
            % Enlarge figure to full screen.
            set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
            % Let's compute and display the histogram.
            [pixelCount, grayLevels] = imhist(grayImage);
            subplot(3, 3, 2); 
            bar(grayLevels, pixelCount);
            grid on;
            title('Histogram of original image', 'FontSize', fontSize);
            xlim([0 grayLevels(end)]); % Scale x axis manually.
            drawnow;
    
    mask = grayImage < 200; % Entire can lid
    % Get rid of holes inside.
    mask = imfill(mask, 'holes');
    % Erode to shrink the mask and get rid of letters.
    se = strel('disk', 45, 0);
    mask = imerode(mask, se);
            subplot(3, 3, 3);
            imshow(mask);
            drawnow;
            title('Mask for can interior', 'FontSize', fontSize);
    % Mask the image
    maskedImage = grayImage;
    maskedImage(~mask) = 0;
            subplot(3, 3, 4);
            imshow(maskedImage);
            title('Masked can lid', 'FontSize', fontSize);
    % Get gradient
    gradImage = imgradient(maskedImage);
            subplot(3, 3, 5);
            imshow(gradImage, []);
            title('Gradient', 'FontSize', fontSize);
    % Shrink mask again to get rid of outer outline
    se = strel('disk', 15, 0);
    mask = imerode(mask, se);
    gradImage(~mask) = 0;
            subplot(3, 3, 6);
            imshow(gradImage, []);
            title('Gradient (smaller)', 'FontSize', fontSize);
            % Get histogram
            pixelCount = hist(gradImage(:), 100);
            % Suppress 0 so we can see the histogram
            pixelCount(1)=0;
            subplot(3, 3, 7); 
            bar(pixelCount);
            grid on;
            title('Histogram of Gradient Image', 'FontSize', fontSize);
            drawnow;
    % Threshold to get a binary image
    binaryImage = gradImage > 45;
    % Get rid of holes inside.
    binaryImage = imfill(binaryImage, 'holes');
            subplot(3, 3, 8);
            imshow(binaryImage, []);
            title('Initial Binary Image', 'FontSize', fontSize);
    % Label the image
    labeledImage = bwlabel(binaryImage);
    % Make measurements of orientation
    measurements = regionprops(labeledImage, 'Orientation', 'Area');
    % Find the largest blob
    [allAreas, sortIndexes] = max(sort([measurements.Area], 'Descend'))
    % Pluck out largest one
    biggestBlob = ismember(labeledImage, sortIndexes(1)) > 0;
            subplot(3, 3, 9);
            imshow(biggestBlob, []);
            title('Final Binary Image', 'FontSize', fontSize);
    % Measure again, this time just the largest blob.
    % Label the image
    labeledImage = bwlabel(biggestBlob);
    % Make measurements of orientation
    measurements = regionprops(labeledImage, 'Orientation', 'Area');
    message = sprintf('The orientation angle = %f degrees\n', measurements(1).Orientation)

end