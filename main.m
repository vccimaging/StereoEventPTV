clear all;
clc;
%% 2D tracking part
addpath('./utils');
data_dir = './data';
data_left_file = 'events_master_votex_2_5_2s.mat';
data_right_file = 'events_slave_votex_2_5_2s.mat';
data_left_file = fullfile(data_dir, data_left_file);
data_right_file = fullfile(data_dir, data_right_file);
cam_resolution = [480, 360];
events_display_time_range = 0.004; %% events time slice
iter_num = 5;
framepeak_left = tracking2d(data_left_file, cam_resolution, events_display_time_range, iter_num);
framepeak_right = tracking2d(data_right_file, cam_resolution, events_display_time_range, iter_num);
%% Output path
result_dir = './result';
if ~exist(result_dir,'dir')
    mkdir(result_dir);
end
result_name = 'uv_result.mat';
result_file = fullfile(result_dir, result_name);
%% frame to process
frame_start = 1;
frame_stop = 3;
frame_num = frame_stop-frame_start+1;
iter_num = frame_num;

%% Load Calibration file

load('./data/resolution.mat');
mappingdistance = 1;
mappingdis=20;
mindis = 40;
tracking2dto3dmapping();

pos_history_left = pos_history_left;
pos_history_right = pos_history_right;
vol_left = vol_left;
vol_right = vol_right;

%% Optimization framework
%%% parameters
maxwarping = 5;
pyramid_levels = 3;
size_xyz_level2 = size_xyz/8;
size_xyz_level3 = size_xyz/64;
[x1,y1,z1] = meshgrid(x_res, y_res, z_res);

%%% build pyramid
s = size(x1);
pyramid_size = zeros( pyramid_levels, 3);
pyramid_size(1, :) = s; %x1 need to be changed
for k = 2 : pyramid_levels
    pyramid_size(k, :) = pyramid_size(k-1,:) .* [0.5 0.5 0.5];
end
uv = zeros([s,3, frame_num]);

smooth_sigma = 0.7;
% smooth_sigma = 1.2;
bound_cond='symmetric';
f3 = fspecial3('gaussian', 2 * round(1.5 * smooth_sigma)+1); % change
position=pyramid_vol{1}.occupancy_left_current;
vol_left_current_map = pyramid_vol{1}.vol_left_current_map;
vol_right_current_map = pyramid_vol{1}.vol_right_current_map;
pyramid{1} = zeros([pyramid_size(1,:), frame_num]);
pyra_vol_left{1} = zeros([pyramid_size(1,:),2,frame_num]);
pyra_vol_right{1} = zeros([pyramid_size(1,:),2,frame_num]);
for i = 1 : frame_num
    pyramid{1}(:,:,:,i) = reshape(position(:,i),s);
    pyramid{1}(:,:,:,i) = imfilter(pyramid{1}(:,:,:,i),f3,bound_cond);
    pyra_vol_left{1}(:,:,:,1,i) = reshape(vol_left_current_map(:,1,i),s);
    pyra_vol_left{1}(:,:,:,1,i) = imfilter(pyra_vol_left{1}(:,:,:,1,i), f3,bound_cond);
    pyra_vol_left{1}(:,:,:,2,i) = reshape(vol_left_current_map(:,2,i),s);
    pyra_vol_left{1}(:,:,:,2,i) = imfilter(pyra_vol_left{1}(:,:,:,2,i), f3,bound_cond);
    pyra_vol_right{1}(:,:,:,1,i) = reshape(vol_right_current_map(:,1,i),s);
    pyra_vol_right{1}(:,:,:,1,i) = imfilter(pyra_vol_right{1}(:,:,:,1,i), f3,bound_cond);
    pyra_vol_right{1}(:,:,:,2,i) = reshape(vol_right_current_map(:,2,i),s);
    pyra_vol_right{1}(:,:,:,2,i) = imfilter(pyra_vol_right{1}(:,:,:,2,i), f3,bound_cond);
end
%%% build pyramid
for k = 2 : pyramid_levels
    pyramid_size(k,:) = pyramid_size(k-1,:) .* [0.5 0.5 0.5];
    pyramid{k} = zeros([pyramid_size(k, :), frame_num]);
    pyra_vol_left{k} = zeros([pyramid_size(k, :), 2, frame_num]);
    pyra_vol_right{k} = zeros([pyramid_size(k, :), 2, frame_num]);
    for i = 1 : frame_num
        pyramid{k-1}(:,:,:,i) = imfilter( pyramid{k-1}(:,:,:,i), f3,bound_cond);
        pyramid{k}(:,:,:,i) = resize( pyramid{k-1}(:,:,:,i) , pyramid_size(k, :));
        
        pyra_vol_left{k-1}(:,:,:,1,i) = imfilter(pyra_vol_left{k-1}(:,:,:,1,i), f3,bound_cond);
        pyra_vol_left{k}(:,:,:,1,i) = resize( pyra_vol_left{k-1}(:,:,:,1,i) , pyramid_size(k, :));
        pyra_vol_left{k-1}(:,:,:,2,i) = imfilter(pyra_vol_left{k-1}(:,:,:,2,i), f3,bound_cond);
        pyra_vol_left{k}(:,:,:,2,i) = resize( pyra_vol_left{k-1}(:,:,:,2,i) , pyramid_size(k, :));
        
        pyra_vol_right{k-1}(:,:,:,1,i) = imfilter(pyra_vol_right{k-1}(:,:,:,1,i), f3,bound_cond);
        pyra_vol_right{k}(:,:,:,1,i) = resize( pyra_vol_right{k-1}(:,:,:,1,i) , pyramid_size(k, :));
        
        pyra_vol_right{k-1}(:,:,:,2,i) = imfilter(pyra_vol_right{k-1}(:,:,:,2,i), f3,bound_cond);
        pyra_vol_right{k}(:,:,:,2,i) = resize( pyra_vol_right{k-1}(:,:,:,2,i) , pyramid_size(k, :));
    end
end
%%% Reconstruction

lambda1 = 2.5e-5;
lambda2 = 0.025;
lambda3 = 2.5e-5;
INNER_ITER = 2;
levels = 3;

[x2,y2,z2] = meshgrid(linspace(x_res(1),x_res(end),size_x/(2^(levels-1))), linspace(y_res(1),y_res(end),size_y/(2^(levels-1))), linspace(z_res(1),z_res(end),size_z/(2^(levels-1))));
position=pyramid{levels};
levelsize=pyramid_size(levels,1)*pyramid_size(levels,2)*pyramid_size(levels,3);
vol_left_current_map = pyra_vol_left{levels};
vol_left_current_map = reshape(vol_left_current_map,[levelsize,2,frame_num]);
vol_right_current_map = pyra_vol_right{levels};
vol_right_current_map = reshape(vol_right_current_map,[levelsize,2,frame_num]);

uv_tmp = zeros([pyramid_size(levels,:),3,frame_num]);

uv = estimatelayerflowdown(position, uv_tmp, lambda1, lambda2, lambda3 ,  ...
    INNER_ITER, maxwarping, x2, y2, z2, vol_left_current_map, vol_right_current_map , trans1, trans2, frame_num, levels,pyramid_levels);


for levels = pyramid_levels-1:-1:1
    [x3,y3,z3] = meshgrid(linspace(x_res(1),x_res(end),size_x/(2^(levels-1))), linspace(y_res(1),y_res(end),size_y/(2^(levels-1))), linspace(z_res(1),z_res(end),size_z/(2^(levels-1))));
    position=pyramid{levels};
    levelsize=pyramid_size(levels,1)*pyramid_size(levels,2)*pyramid_size(levels,3);
    vol_left_current_map = pyra_vol_left{levels};
    vol_left_current_map = reshape(vol_left_current_map,[levelsize,2,frame_num]);
    vol_right_current_map = pyra_vol_right{levels};
    vol_right_current_map = reshape(vol_right_current_map,[levelsize,2,frame_num]);
    
    uv_tmp = zeros([pyramid_size(levels,:),3,frame_num]);
    for i = 1 : frame_num
        uv_tmp(:,:,:,:,i) = cat(4,interp_new(x2,y2,z2,uv(:,:,:,1,i),x3,y3,z3,'cubic'),...
            interp_new(x2,y2,z2,uv(:,:,:,2,i),x3,y3,z3,'cubic'),...
            interp_new(x2,y2,z2,uv(:,:,:,3,i),x3,y3,z3,'cubic'));
    end
    uv = estimatelayerflowdown(position,  uv_tmp, lambda1,  lambda2, lambda3 ,...
        INNER_ITER, maxwarping, x3, y3, z3, vol_left_current_map, vol_right_current_map , trans1, trans2, frame_num, levels,pyramid_levels);
    x2=x3;
    y2=y3;
    z2=z3;
end

%% Result save


save(result_file,'uv','x_res','y_res','z_res','-v7.3')

