function test_detect_noise()

    noise_sigma = 10;
    noise_hsize = [10 10];
    noise_threshold = 1;

    noise_im = zloadim('movies/Timelapse 2-02(9)_mostlyunsat_s9.ome-NOISE_TEST.tif');
%     signal_im = zloadim('movies/single_molecule_signal.tif');
    
    [noise_im_is_noise,noise_im_filtered_list,noise_im_filtered] = detect_noise( noise_im, noise_threshold, noise_hsize, noise_sigma );
%     [signal_im_is_noise,signal_im_filtered_list,signal_im_filtered] = detect_noise( signal_im, noise_threshold, noise_hsize, noise_sigma );
    
    for frame_idx = 1:size(noise_im,3)

        figure(130892)
        clf
        
        % Original images
        
            subplot(3,2,1)
                hold all

                imagesc(noise_im(:,:,frame_idx))
                
                colormap gray
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                do_box(noise_im_is_noise(frame_idx))
                
                title('(Nearly) pure noise image')

%             subplot(3,2,2)
%                 hold all
% 
%                 imagesc(signal_im(:,:,frame_idx))
%                 
%                 colormap gray
%                 axis image
%                 set(gca,'Ydir','Reverse')
%                 axis off
%                 
%                 do_box(signal_im_is_noise(frame_idx))
%                 
%                 title('Signal image')

        % Entropy images
            cur_noise_filtered = noise_im_filtered(:,:,frame_idx);
%             cur_signal_filtered = signal_im_filtered(:,:,frame_idx);
                
            subplot(3,2,3)
                hold all

                imagesc(cur_noise_filtered)

                colormap gray
                axis image
                set(gca,'Ydir','Reverse')
                axis off
                
                do_box(noise_im_is_noise(frame_idx))
                
                title('Entropy-filtered noise image')

%             subplot(3,2,4)
%                 hold all
% 
%                 imagesc(cur_signal_filtered)
%                 
%                 colormap gray
%                 axis image
%                 set(gca,'Ydir','Reverse')
%                 axis off
%                 
%                 do_box(signal_im_is_noise(frame_idx))
% 
%                 title('Entropy-filtered signal image')
                
        % Entropy histograms
        
            subplot(3,2,5)
                hold all
                
                nonzero = cur_noise_filtered(cur_noise_filtered ~= 0);
                hist(nonzero(:),70)
                
                title('Noise entropy distribution')
                xlabel('Normalized entropy')
                
%             subplot(3,2,6)
%                 hold all
%                 
%                 nonzero = cur_signal_filtered(cur_signal_filtered ~= 0);
%                 hist(nonzero(:),70)
%                 
%                 xl = xlim;
%                 yl = ylim;
%                 
%                 title('Signal entropy distribution')
%                 xlabel('Normalized entropy')
                
%             subplot(3,2,5)
%                 xlim(xl)
%                 ylim(yl)

        set(findall(gcf,'type','text'),'fontSize',20,'fontWeight','bold')
        set(findall(gcf,'type','axes'),'fontSize',20,'fontWeight','bold','LineWidth',3)
        set(gcf, 'color', 'white');

        drawnow
    end
end

function do_box(is_noise)
    xl = xlim;
    yl = ylim;
        
    if is_noise
        line(xlim,[yl(1) yl(1)],'LineWidth',3,'Color','r')
        line(xlim,[yl(2) yl(2)],'LineWidth',3,'Color','r')
        line([xl(1) xl(1)], ylim,'LineWidth',3,'Color','r')
        line([xl(2) xl(2)], ylim,'LineWidth',3,'Color','r')
    else
        line(xlim,[yl(1) yl(1)],'LineWidth',3,'Color','g')
        line(xlim,[yl(2) yl(2)],'LineWidth',3,'Color','g')
        line([xl(1) xl(1)], ylim,'LineWidth',3,'Color','g')
        line([xl(2) xl(2)], ylim,'LineWidth',3,'Color','g')
    end
end