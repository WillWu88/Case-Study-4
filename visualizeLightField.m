%% Case study 4: Ray tracing in optics
% *ESE 105* 
%
% *Name: FILL IN HERE*

clear;
close all;

% load light field
load('lightField.mat');

% visualize rays as an image
[rayImg,x,y] = rays2img(ray_x,ray_y,ray_color,.025,300);
image(x,y,rayImg); axis image xy;
xlabel('x (m)'); ylabel('y (m)');
