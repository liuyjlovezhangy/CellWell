function [ output_args ] = test_tm( input_args )
%TEST_TM Summary of this function goes here
%   Detailed explanation goes here

    % Find maximum response
    I = im2double(imread('big_img.tif'));
    % Template of Eye Lena
    T = im2double(imread('well4.tif'));
    % Calculate SSD and NCC between Template and Image
    [I_SSD,I_NCC]=template_matching(T,I);

    num_levels = 13;
    
    [IDX,sep] = otsu(I_NCC,num_levels);
    sep
    
    max_map = IDX ==num_levels;
    
%     IDX==num_levels
    
%     find(IDX==num_levels)
    
    figure(15432)
    clf

    subplot(1,2,1)
        hold all
        hist(I_NCC(:),2000)
    %    hist(I_SSD(:),2000)
    
    subplot(1,2,2)
        hold all
        
%         imagesc(IDX)
        imagesc(max_map)
        
        colormap gray
        axis tight
        axis equal
        
    return

    % Find maximum correspondence in I_SDD image
    [x,y]=find(I_SSD==max(I_SSD(:)));
    % Show result
    figure, 
    subplot(2,2,1), imshow(I); hold on; plot(y,x,'r*'); title('Result')
    subplot(2,2,2), imshow(T); title('The eye template');
    subplot(2,2,3), imshow(I_SSD); title('SSD Matching');
    subplot(2,2,4), imshow(I_NCC); title('Normalized-CC');

return

  % Find maximum response
   I = im2double(imread('lena.jpg'));
  % Template of Eye Lena
   T=I(124:140,124:140,:);
  % Calculate SSD and NCC between Template and Image
   [I_SSD,I_NCC]=template_matching(T,I);
  % Find maximum correspondence in I_SDD image
   [x,y]=find(I_SSD==max(I_SSD(:)));
  % Show result
   figure, 
   subplot(2,2,1), imshow(I); hold on; plot(y,x,'r*'); title('Result')
   subplot(2,2,2), imshow(T); title('The eye template');
   subplot(2,2,3), imshow(I_SSD); title('SSD Matching');
   subplot(2,2,4), imshow(I_NCC); title('Normalized-CC');
end

