function frameinfo = tracking2d(filename, cam_resolution, events_display_time_range, iter_num)
defineparams();
% load(filename);
curr_events = load(filename);
fields = fieldnames(curr_events);
curr_events = curr_events.(fields{1});
curr_events(:,4) = 1;
curr_events_pos = curr_events;
curr_events_neg = curr_events;
% curr_events_pos(:,4)=1;

curr_events_pos = curr_events(curr_events(:,4) > 0, :);
curr_events_neg = curr_events(curr_events(:,4) < 0, :);
curr_events_pos_mask = zeros(size(curr_events_pos, 1), 1);
curr_events_neg_mask = zeros(size(curr_events_neg, 1), 1);

rows = cam_resolution(1);
cols = cam_resolution(2);


event_t0 = 0;
last_event_iter = 1;
start_event_iter = find(event_t0 <= curr_events_pos(:, 1), 1);
new_event_iter = find(event_t0 + events_display_time_range <= curr_events_pos(:,1),1);
old_event_iter = new_event_iter;
curr_event_display = curr_events_pos(start_event_iter : new_event_iter -1, :);

event_t0 = curr_events_pos(new_event_iter - 1,1);

neg_all = curr_event_display(:,4) < 0;
curr_event_display_pos = curr_event_display( curr_event_display(:,4) > 0, :);
curr_event_display_neg = curr_event_display( curr_event_display(:,4) < 0, :);
curr_events_pos_image{1} = accumarray(round([curr_event_display_pos(:,2) curr_event_display_pos(:, 3)]) + 1, 1, [rows cols]); %% accumulated images


%% Particle Location Identification
filter_size = 6;
img_window = 8;
filter = zeros(filter_size,filter_size);
for i = 1 : filter_size
    for j = 1 : filter_size
        if sqrt((i -0.5 - filter_size/2)^2 + (j - 0.5- filter_size/2)^2) <= 3
            filter(i,j) = 1;
        end
    end
end
peak_find_size2 = [12,12];
curr_events_pos_image2{1} = curr_events_pos_image{1};
curr_events_pos_image2{1}(curr_events_pos_image2{1} >0) = 1;
img2 = imfilter(curr_events_pos_image{1}, filter);
%% Particle center detection
peak2 = pkfnd(img2, peak_find_size2(1), peak_find_size2(2));
peak_center = zeros(size(peak2));
for i = 1 : size(peak2,1)
    peak_window_index = img2(peak2(i,1) - img_window/2 : peak2(i,1) + img_window/2, ...
        peak2(i,2)-img_window/2: peak2(i,2) + img_window/2);
    peak_win_sum = sum(peak_window_index(:));
    [y_mesh,x_mesh] = meshgrid( peak2(i,2) - img_window/2:peak2(i,2) + img_window/2,peak2(i,1) - img_window/2:peak2(i,1) + img_window/2);
    location_x = peak_window_index .* x_mesh/peak_win_sum;
    location_y = peak_window_index .* y_mesh/peak_win_sum;
    peak_center(i,:) = [sum(location_x(:)), sum(location_y(:))];
    
end
peak2 = peak_center;

%% particle center checking
peak_mark = false(length(peak2),1);
events_marker = false(length(curr_event_display_pos),1);
vol_ini = zeros(size(peak2));
flow_ini = [0;0];
for i = 1 : length(peak2)
    %     i  = 194
    dt = curr_event_display_pos(end, 1) - curr_event_display_pos(1,1);
    [flow, shifted_points, event_marker] = em_flow(...
        curr_event_display_pos,...
        peak2(i,:),...
        flow_ini,...
        params);
    events_marker(event_marker) = true;
    if sum(isnan(flow)) > 0
        peak_mark(i) = false;
    else
        peak_mark(i) = true;
        vol_ini(i,:) = flow'*dt;
    end
end
peak_tmp = peak2(peak_mark,:);
frameinfo(1).peak=peak_tmp;
frameinfo(1).time = dt;
frameinfo(1).vol = vol_ini(peak_mark,:);
vol_before = vol_ini(peak_mark,:);
peaks_before = peak_tmp;
peaks_predict = peak_tmp + vol_ini(peak_mark,:);

curr_events_pos_mask(start_event_iter : new_event_iter - 1) = events_marker;
iter_num_max = floor(curr_events_pos(end,1)/events_display_time_range) - 1;
if (iter_num > iter_num_max)
    fprintf('Frame number exist the frame number in the dataset');
    return;
end
fprintf('2D tracking frame: ');
for iter = 1 : iter_num
    fprintf('%d..',iter);
    int_time = events_display_time_range;
    new_event_iter = find(int_time + event_t0 <= curr_events_pos(:,1),1);
    curr_event_display = curr_events_pos(old_event_iter : new_event_iter - 1,:);
    event_t0_old = event_t0;
    event_t0 = curr_events_pos(new_event_iter -1,1);
    curr_event_display_pos = curr_event_display;
    curr_events_pos_image{iter + 1} = accumarray(round([curr_event_display_pos(:,2) curr_event_display_pos(:, 3)]) + 1, 1, [rows cols]);
    images_test = curr_events_pos_image{iter + 1};
    images_test(images_test>1) =1;
    img2 = imfilter(images_test, filter);
    peak2 = pkfnd(img2, peak_find_size2(1), peak_find_size2(2));
    peak_center = zeros(size(peak2));
    for i = 1 : size(peak2,1)
        peak_window_index = img2( peak2(i,1) - img_window/2 : peak2(i,1) + img_window/2, peak2(i,2) - img_window/2 : peak2(i,2) + img_window/2);
        peak_win_sum = sum(peak_window_index(:));
        [y_mesh,x_mesh] = meshgrid( peak2(i,2) - img_window/2:peak2(i,2) + img_window/2,peak2(i,1) - img_window/2:peak2(i,1) + img_window/2);
        location_x = peak_window_index .* x_mesh/peak_win_sum;
        location_y = peak_window_index .* y_mesh/peak_win_sum;
        peak_center(i,:) = [sum(location_x(:)), sum(location_y(:))];
    end
    
    peak2 = peak_center;
    peak_mark = false(length(peak2),1);
    events_marker = false(length(curr_event_display_pos),1);
    vol_ini = zeros(size(peak2));
    flow_ini = [0;0];
    for i = 1 : length(peak2)
        dt = curr_event_display_pos(end,1) - curr_event_display_pos(1,1);
        [flow, shifted_points, event_marker] = em_flow(...
            curr_event_display_pos,...
            peak2(i,:),...
            flow_ini,...
            params);
        events_marker(event_marker) = true;
        if sum(isnan(flow)) > 0
            peak_mark(i) = false;
        else
            peak_mark(i) = true;
            vol_ini(i,:) = flow'*dt;
        end
    end
    curr_events_pos_mask(old_event_iter : new_event_iter -1) = events_marker;
    old_event_iter = new_event_iter;
    peak_tmp = peak2(peak_mark, :);
    vol_tmp = vol_ini(peak_mark, :);
    frameinfo(iter+1).peak = peak_tmp;
    frameinfo(iter+1).time = dt;
    frameinfo(iter+1).vol = vol_tmp;
    peak_predict_back = peak2(peak_mark,:) - vol_ini(peak_mark,:);
    
    peak_dist = (peaks_predict(:,1) - peak_tmp(:,1)').^2 + (peaks_predict(:,2) - peak_tmp(:,2)').^2;
    peak_dist2 = (peak_predict_back(:,1)' - peaks_before(:,1)).^2 + (peak_predict_back(:,2)' - peaks_before(:,2)).^2;
    
    vol_before_max = 8*(vol_before(:,1).^2 + vol_before(:,2).^2);
    vol_tmp_max = 8*(vol_tmp(:,1).^2 + vol_tmp(:,2).^2);
    
    for kk = 1 : size(vol_before,1)
        peak_dist( kk,peak_dist(kk,:) > vol_before_max(kk)) = inf;
    end
    for kk = 1 : size(vol_tmp, 1)
        peak_dist2( peak_dist2(:,kk) > vol_tmp_max(kk), kk) = inf;
    end
    
    
    %     peak_dist(peak_dist>60) = inf;
    [assignment,cost] = munkres(peak_dist);
    
    %     peak_dist2(peak_dist2>60) = inf;
    [assignment2,cost2] = munkres(peak_dist2);
    
    assignment_final = assignment;
    assignment_index = assignment_final == 0;
    assignment_final(assignment_index)=assignment2(assignment_index);
    
    peak_next = zeros(size(peaks_predict,1),2);
    index = assignment_final~=0;
    ind = 1 : length(assignment_final);
    indx = ind(index);
    
    peak_next(index,:) = peak_tmp(assignment_final(index),:);
    vel_rel = peak_tmp(assignment_final(index),:) -peaks_before(index,:);
    vel_comp = vol_before(index,:);
    vel_cor = sum(vel_rel.*vel_comp,2);
    index(indx(vel_cor<0)) = false;  %% false detection
    
    vel_rel = peak_tmp(assignment_final(index),:) -peaks_before(index,:);
    frameinfo(iter).peak_val = peaks_before(index,:);
    frameinfo(iter).vel_est = [vel_rel; vol_before(~index,:)];
    %     frameinfo(events_iteration-1).vel_est_modi = vel_rel2;
    
    frameinfo(iter).volmap = vol_before(index,:);
    frameinfo(iter).peak_current = [peaks_before(index,:); peaks_before(~index,:)];
    frameinfo(iter).peak_next=peak_next(index,:);
    peaks_predict = peak_tmp + vol_ini(peak_mark,:);
    peaks_before = peak_tmp;
    vol_before = vol_ini(peak_mark,:);
    
    
end
fprintf('\n');
end




