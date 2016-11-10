function improved_rotation()
    angle_manual = -2.5;

    se2 = strel('line',10,90);
    se3 = strel('disk',50);
    
    im = imread('C5-Microwells_4color-04(10)_s10_BRIGHT-1.tif');
    im = mat2gray(im);
    
    im = wiener2(im);
    
    im2 = im;
    
    
%     im2 = imadjust(im2,[],[],0.50);
%     im2 = imresize(im2,2);

    im2 = imopen(im2,se2);
    
    im3 = imtophat(im,se3);
    
    angle1 = horizon(im,0.1,'hough')
    angle2 = horizon(im2,0.1,'hough')
    angle3 = horizon(im3,0.1,'hough')
    
    im_rotato = imrotate(im,-angle1,'bicubic','crop');
    im2_rotato = imrotate(im2,-angle2,'bicubic','crop');
    im3_rotato = imrotate(im3,-angle3,'bicubic','crop');
    
    im_rotato_width = size(im_rotato,2);
    im_rotato_height = size(im_rotato,1);
    
    figure(1432)
    clf
    
        subtightplot(2,3,1)
            hold all
        
            imagesc(im)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 

            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')
            
            colormap gray
            
            freezeColors
            
        subtightplot(2,3,2)
            hold all
        
            imagesc(im2)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 
            
            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,3)
            hold all
        
            imagesc(im3)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 
            
            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,4)
            hold all
        
            imagesc(im_rotato)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 
            
            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')

            colormap gray
            
            freezeColors
        
        subtightplot(2,3,5)
            hold all
        
            imagesc(im2_rotato)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 
            
            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')

            colormap gray
            
            freezeColors
            
        subtightplot(2,3,6)
            hold all
        
            imagesc(im3_rotato)
            
            axis image
            set(gca,'Ydir','Reverse')
            axis off 
            
            line(im_rotato_width/2*[1 1],ylim,'LineStyle','--','LineWidth',5,'Color','r')
            line(xlim,im_rotato_height/2*[1 1],'LineStyle','--','LineWidth',5,'Color','r')

            colormap gray
            
            freezeColors
end

