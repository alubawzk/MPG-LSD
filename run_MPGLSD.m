clc;clear;
close all;
%%
% mex mpglsd.c

%% image root
Noise_Level = 'Reference';
% Noise_Level = 'Noisy';

% York
% load('D:/360MoveData/Users/ylg/Desktop/LineDetection/DataSets/YorkUrban-LineSegment/LineSegmentAnnotation/Image_ID_List.mat');
% img_root = ['D:/360MoveData/Users/ylg/Desktop/LineDetection/DataSets/YorkUrban-LineSegment/',Noise_Level];
% output_root = ['D:/360MoveData/Users/ylg/Desktop/LineDetection/Evaluations/InputData/MSLSD/' ,Noise_Level, '/'];
% NoI = 102; type = '.jpg';

% Renior
% load('D:/360MoveData/Users/ylg/Desktop/LineDetection/DataSets/Renoir-LineSegment/Reference/LineSegmentAnnotation/Image_ID_List.mat');
% img_root = ['D:/360MoveData/Users/ylg/Desktop/LineDetection/DataSets/Renoir-LineSegment/',Noise_Level];
% output_root = ['D:/360MoveData/Users/ylg/Desktop/LineDetection/Evaluations_Bigscale/InputData/MSLSD/' ,Noise_Level, '/'];
% NoI = 40; type = '.bmp';

%%
flag = 1; % 0: dataset, 1: single image
scale_num = 5;
condition = 3; %3

%% run
if flag % single image
%     i_im = 11;
%     img_name = Image_ID_List(i_im).name;
%     img_path = sprintf('%s/%s/%s.jpg', img_root, img_name, img_name);
%     input_img  = imread(img_path);
    input_img  = imread('P1020829.jpg');
    im = input_img;
    im_gray = rgb2gray(im);
    [lines] = MPGLSD( im_gray, 1, scale_num, condition );
else % dataset
    for i_im = 1:NoI
        if i_im<10
            fprintf('   %d',i_im)
        elseif i_im<100
            fprintf('  %d',i_im)
        else
            fprintf(' %d',i_im)
        end
        if mod(i_im,17)==0
            fprintf('\n')
        end
        img_name = Image_ID_List(i_im).name;
        img_path = sprintf(['%s/%s/%s' type], img_root, img_name, img_name);
        input_img  = imread(img_path);
        im = input_img;
        im_gray = rgb2gray(im);
        [lines] = MPGLSD( im_gray, 0, scale_num, condition );
        output_dir = [output_root, 'im' num2str(i_im)];
        if exist(output_dir,'dir')==0
            mkdir(output_dir);
        end
        lineset = lines(:,1:5);
        save([output_dir,'/literature.mat'], 'lineset');
    end
end
