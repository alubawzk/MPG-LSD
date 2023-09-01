clc;clear;
close all;

%% run
draw_flag  = 1;
input_img  = imread('P1020829.jpg');
im = input_img;
im_gray = rgb2gray(im);
[lines] = MPGLSD( im_gray, draw_flag);
