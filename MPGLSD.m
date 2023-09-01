function lines = MPGLSD(img, draw_flag)  
[lines,~,~,~] = mpglsd(double(img));
lines = lines(:,1:5);

%% »­Í¼
if draw_flag
    [m,n] = size(img(:,:,1));
    ls = lines(:,1:4) ;
    figure;
    imshow(ones(m,n)*255,'border','tight');
    hold on;
    for i = 1:length(ls(:,1))
        plot([ls(i,1) ls(i,3)],[ls(i,2) ls(i,4)], 'k','linewidth',0.75 );
    end
    plot([1 n n 1 1],[1 1 m m 1],'color',[.75 .75 .75],'linewidth',0.75);
    hold off
end  

%% draw gt
% load('LineSegmentAnnotation/Image_ID_List.mat');
% str_gnd = sprintf('LineSegmentAnnotation/%s_GND.mat', Image_ID_List(i_im).name);
% load(str_gnd);
% ls_gnd = unique(line_gnd, 'rows');
% if draw_flag
%     [m,n] = size(img(:,:,1));
% %     figure;
%     subplot(1,2,2)
%     imshow(ones(m,n)*255,'border','tight');
%     hold on;
%     for i = 1:length(ls_gnd(:,1))
%         plot([ls_gnd(i,1) ls_gnd(i,3)],[ls_gnd(i,2) ls_gnd(i,4)], 'k','linewidth',1 );
%     end
%     plot([1 n n 1 1],[1 1 m m 1],'color',[.75 .75 .75],'linewidth',0.75);
%     hold off
% end 
end

