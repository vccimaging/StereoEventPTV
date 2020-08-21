%% Extract data from 2D tracking
peaknum = [];
for i = 1 : iter_num
    peaknum = [peaknum;[length(framepeak_left(i).peak),length(framepeak_right(i).peak)]];
end
traj_num = max(peaknum(:));
pos_history_left = zeros(traj_num,2,iter_num);
vol_left = zeros(traj_num,2,iter_num);
pos_history_right = zeros(traj_num,2,iter_num);
vol_right = zeros(traj_num,2,iter_num);

for i = 1 : iter_num
    pos_history_left(1:peaknum(i,1),:,i) = framepeak_left(i).peak_current;
    vol_left(1:peaknum(i,1),:,i) = framepeak_left(i).vel_est;
    pos_history_right(1:peaknum(i,2),:,i) = framepeak_right(i).peak_current;
    vol_right(1:peaknum(i,2),:,i) = framepeak_right(i).vel_est;
end


T3_matrix = [ trans1(1,1), trans1(2,1), trans1(3,1);...
    trans1(1,2), trans1(2,2), trans1(3,2);...
    trans2(1,1), trans2(2,1), trans2(3,1);...
    trans2(1,2), trans2(2,2), trans2(3,2)];

T4_matrix = inv(T3_matrix' * T3_matrix) *T3_matrix';


size_z = length(z_res);

tol = 4;
%% Calculating voxel boudaries
index_z = 1 : size_z;
eval_points = [0,0;480,0;0,360;480,360];
boundaries = zeros(size_z,4);
for i = 1 : size_z
    point_val = zeros(2,8);
    for j = 1 : 4
        A_matrix = [trans1(1,1)-trans1(1,3)*eval_points(j,1), trans1(2,1) - trans1(2,3)*eval_points(j,1);...
            trans1(1,2)-trans1(1,3)*eval_points(j,2), trans1(2,2) - trans1(2,3)*eval_points(j,2)];
        b_matrix = [trans1(3,3)*eval_points(j,1)*z_res(i) + eval_points(j,1) - trans1(3,1) * z_res(i) - trans1(4,1);...
            trans1(3,3)*eval_points(j,2)*z_res(i) + eval_points(j,2) - trans1(3,2) * z_res(i) - trans1(4,2)];
        point_mapping1 = inv(A_matrix) * b_matrix;
        A_matrix = [trans2(1,1)-trans2(1,3)*eval_points(j,1), trans2(2,1) - trans2(2,3)*eval_points(j,1);...
            trans2(1,2)-trans2(1,3)*eval_points(j,2), trans2(2,2) - trans2(2,3)*eval_points(j,2)];
        b_matrix = [trans2(3,3)*eval_points(j,1)*z_res(i) + eval_points(j,1) - trans2(3,1) * z_res(i) - trans2(4,1);...
            trans2(3,3)*eval_points(j,2)*z_res(i) + eval_points(j,2) - trans2(3,2) * z_res(i) - trans2(4,2)];
        point_mapping2 = inv(A_matrix) * b_matrix;
        point_val(:,2*(j-1)+1) = point_mapping1';
        point_val(:,2*j) = point_mapping2';
    end
    point_val = point_val + 10;
    boundaries(i,1) = max([point_val(1,1),point_val(1,2), point_val(1,5),point_val(1,6)]) -10;
    boundaries(i,2) = max([point_val(2,1),point_val(2,2), point_val(2,3),point_val(2,4)]) -10;
    boundaries(i,3) = min([point_val(1,3),point_val(1,4), point_val(1,7),point_val(1,8)]) -10;
    boundaries(i,4) = min([point_val(2,5),point_val(2,6), point_val(2,7),point_val(2,8)]) -10;
end



mindis = 5;
parapairs=cell(frame_num,1);
indexvalue = 1:size(pos_history_left,1);
%% 3D mapping
fprintf('processing frame ');
for frame = frame_start : frame_stop
    value=0.1;
    % frame = 1;
    fprintf('%d..',frame);
    frame_idx = frame - frame_start+1;
    pos_left_current = pos_history_left(:, :,frame_idx);
    vol_left_current = vol_left(:,:,frame_idx);
    
    pos_right_current = pos_history_right(:,:,frame_idx);
    vol_right_current = vol_right(:,:,frame_idx);
    
    pos_marker = pos_left_current(:,1)~=0;
    pos_left_current = pos_left_current(pos_marker,:);
    vol_left_current = vol_left_current(pos_marker,:);
    
    pos_marker = pos_right_current(:,1) ~=0;
    pos_right_current = pos_right_current(pos_marker, :);
    vol_right_current = vol_right_current(pos_marker, :);
    
    pos_left_size = size(pos_left_current,1);
    pos_right_size = size(pos_right_current,1);
    index_left = 1 : pos_left_size;
    index_right = 1 : pos_right_size;
    occup_matrix_left = zeros(pos_left_size,pos_right_size);
    occup_matrix_right = zeros(pos_left_size,pos_right_size);
    dist_matrix_left = zeros(pos_left_size,pos_right_size);
    dist_matrix_right= zeros(pos_left_size,pos_right_size);
    min_dis_left = [];
    min_dis_right = [];
    for i = 1 : pos_left_size
        A_matrix = [trans1(1,1)-trans1(1,3)*pos_left_current(i,1), trans1(2,1) - trans1(2,3)*pos_left_current(i,1);...
            trans1(1,2)-trans1(1,3)*pos_left_current(i,2), trans1(2,2) - trans1(2,3)*pos_left_current(i,2)];
        b_matrix = [trans1(3,3)*pos_left_current(i,1)*z_res + pos_left_current(i,1) - trans1(3,1) * z_res - trans1(4,1);...
            trans1(3,3)*pos_left_current(i,2)*z_res + pos_left_current(i,2) - trans1(3,2) * z_res - trans1(4,2)];
        point_mapping = inv(A_matrix) * b_matrix;
        point_mapping_index = point_mapping(1,:) >= boundaries(:,1)' & point_mapping(1,:) <= boundaries(:,3)' & ...
            point_mapping(2,:) >= boundaries(:,2)' & point_mapping(2,:) <= boundaries(:,4)';
        point_inbox = point_mapping(:,point_mapping_index);
        
        point_3d = [point_inbox;z_res(point_mapping_index);ones(1,sum(point_mapping_index))]';
        point_right=point_3d * trans2;
        point_right = point_right(:,1:2)./point_right(:,3);
        if size(point_right,1) > 2
            pf = polyfit(point_right(:,1),point_right(:,2),1);
            vec2 = [pos_right_current(:,1) - point_right(1,1), pos_right_current(:,2) - point_right(1,2)];
            vec1 = [point_right(end,1) - point_right(1,1), point_right(end,2) - point_right(1,2)];
            d = vec1(1)^2 + vec1(2)^2;
            min_dist = zeros(size(vec2, 1),1);
            r = (vec1(1)*vec2(:,1) + vec1(2)*vec2(:,2))/d;
            for j = 1 : size(vec2, 1)
                if r(j) <= 0
                    min_dist(j) = sqrt((pos_right_current(j,1) - point_right(1,1))^2 + (pos_right_current(j,2) - point_right(1,2))^2);
                elseif r(j) >=1
                    min_dist(j) = sqrt((pos_right_current(j,1) - point_right(end,1))^2 + (pos_right_current(j,2) - point_right(end,2))^2);
                else
                    theta = (vec2(j,:) * vec1')/(norm(vec2(j,:)) * norm(vec1));
                    if theta >=1
                        theta = 1;
                    elseif theta <=-1
                        theta = -1;
                    end
                    
                    min_dist(j) = sqrt((pos_right_current(j,1) - point_right(1,1))^2 + (pos_right_current(j,2) - point_right(1,2))^2) * sin(acos(theta));
                end
            end
            index1 = min_dist < mindis;
            occup_matrix_left(i, index1) = 1;
            dist_matrix_left(i, index1) = min_dist(index1);
            dist_matrix_left(i, ~index1) = inf;
            
        end
        
    end
    for i = 1 : pos_right_size
        A_matrix = [trans2(1,1) - trans2(1,3)*pos_right_current(i,1), trans2(2,1)-trans2(2,3)*pos_right_current(i,1);...
            trans2(1,2) - trans2(1,3)*pos_right_current(i,2), trans2(2,2)-trans2(2,3)*pos_right_current(i,2)];
        b_matrix = [trans2(3,3)*pos_right_current(i,1)*z_res + pos_right_current(i,1) - trans2(3,1)*z_res - trans2(4,1);...
            trans2(3,3)*pos_right_current(i,2)*z_res + pos_right_current(i,2) - trans2(3,2)*z_res - trans2(4,2)];
        point_mapping = inv(A_matrix) * b_matrix;
        point_mapping_index = point_mapping(1,:) >= boundaries(:,1)' & point_mapping(1,:) <= boundaries(:,3)' & ...
            point_mapping(2,:) >= boundaries(:,2)' & point_mapping(2,:) <= boundaries(:,4)';
        point_inbox = point_mapping(:,point_mapping_index);
        
        point_3d = [point_inbox;z_res(point_mapping_index);ones(1,sum(point_mapping_index))]';
        
        point_left = point_3d * trans1;
        point_left = point_left(:,1:2)./point_left(:,3);
        if size(point_left, 1) > 2
            pf = polyfit(point_left(:,1),point_left(:,2),1);
            vec2 = [pos_left_current(:,1) - point_left(1,1), pos_left_current(:,2) - point_left(1,2)];
            vec1 = [point_left(end,1) - point_left(1,1), point_left(end,2) - point_left(1,2)];
            d = vec1(1)^2 + vec1(2)^2;
            min_dist = zeros(size(vec2,1),1);
            r = (vec1(1)*vec2(:,1) + vec1(2)*vec2(:,2))/d;
            for j = 1 : size(vec2,1)
                if r(j) <= 0
                    min_dist(j) = sqrt((pos_left_current(j,1) - point_left(1,1))^2 + (pos_left_current(j,2) - point_left(1,2))^2);
                elseif r(j) >=1
                    min_dist(j) = sqrt((pos_left_current(j,1) - point_left(end,1))^2 + (pos_left_current(j,2) - point_left(end,2))^2);
                else
                    theta = vec2(j,:) * vec1'/(norm(vec2(j,:)) * norm(vec1));
                    if theta >=1
                        theta = 1;
                    elseif theta <=-1
                        theta = -1;
                    end
                    min_dist(j) = sqrt((pos_left_current(j,1) - point_left(1,1))^2 + (pos_left_current(j,2) - point_left(1,2))^2) * sin(acos(theta));
                end
            end
            
            index2 = min_dist < mindis;
            occup_matrix_right(index2 ,i) = 1;
            dist_matrix_right(index2 ,i) = min_dist(index2 );
            dist_matrix_right(~index2, i) = inf;
        end
        
    end
    occup = occup_matrix_left & occup_matrix_right;
    dist_matrix_left(~occup) = inf;
    dist_matrix_right(~occup) = inf;
    dist_total = dist_matrix_left + dist_matrix_right;
    
    dist_hugs = dist_total;
    
    index_rank = [1,1];
    while(~isempty(index_rank))
        [assignment, cost] = munkres(dist_hugs');
        pairs = [];
        index_agw = [];
        for i = 1 : pos_right_size
            if assignment(i) ~=0
                index_agw = [index_agw, i];
                A_matrix = [trans2(1,1) - trans2(1,3)*pos_right_current(i,1), trans2(2,1)-trans2(2,3)*pos_right_current(i,1);...
                    trans2(1,2) - trans2(1,3)*pos_right_current(i,2), trans2(2,2)-trans2(2,3)*pos_right_current(i,2)];
                b_matrix = [trans2(3,3)*pos_right_current(i,1)*z_res + pos_right_current(i,1) - trans2(3,1)*z_res - trans2(4,1);...
                    trans2(3,3)*pos_right_current(i,2)*z_res + pos_right_current(i,2) - trans2(3,2)*z_res - trans2(4,2)];
                point_mapping = inv(A_matrix) * b_matrix;
                point_mapping_index = point_mapping(1,:) >= boundaries(:,1)' & point_mapping(1,:) <= boundaries(:,3)' & ...
                    point_mapping(2,:) >= boundaries(:,2)' & point_mapping(2,:) <= boundaries(:,4)';
                
                idx = index_z(point_mapping_index);
                point_inbox = point_mapping(:,point_mapping_index);
                point_3d = [point_inbox;z_res(point_mapping_index);ones(1,sum(point_mapping_index))]';
                point_3d_vec = point_3d(end,1:3) - point_3d(1,1:3);
                point_left = point_3d * trans1;
                point_left = point_left(:,1:2)./point_left(:,3);
                vec2 = [pos_left_current(assignment(i),1) - point_left(1,1), pos_left_current(assignment(i),2) - point_left(1,2)];
                vec1 = [point_left(end,1) - point_left(1,1), point_left(end,2) - point_left(1,2)];
                d = vec1(1)^2 + vec1(2)^2;
                r = (vec1(1)*vec2(1) + vec1(2) * vec2(2))./d;
                r1 = round(r * (size(point_3d,1) -1))/(size(point_3d,1)-1);
                point_3d_cos_current = point_3d(1,1:3) + r * point_3d_vec;
                
                
                vel = (T4_matrix * [vol_left_current(assignment(i),:)';vol_right_current(i,:)'])';
                pairs=[pairs;pos_left_current(assignment(i),:),pos_right_current(i,:),point_3d_cos_current,vel,dist_hugs( assignment(i),i), vol_left_current(assignment(i),:), vol_right_current(i,:)];
            end
        end
        pair_points_3d = pairs(:,5:7);
        point_mapping_index = point_mapping(1,:) >= boundaries(:,1)' & point_mapping(1,:) <= boundaries(:,3)' & ...
            point_mapping(2,:) >= boundaries(:,2)' & point_mapping(2,:) <= boundaries(:,4)';
        vel_pair = pairs(:,8:10);
        vel_strength = vecnorm(vel_pair,2,2);
        mean_vel_strength = mean(vel_strength);
        
        %
        % % %   AGW %% To remove critical mismatch
        particle_distance = sqrt((pairs(:,5) - pairs(:,5)').^2+(pairs(:,6) - pairs(:,6)').^2+(pairs(:,7) - pairs(:,7)').^2);
        total_dis = 0;
        pair_size = size(pairs,1);
        neighbor_num = 8;
        delta_w = 0.3;
        lambda_weight = zeros(1,neighbor_num);
        %     index_agw = 1:pair_size;
        index_pair = 1 : pair_size;
        flag_agw = true(pair_size, 1);
        rvalue = zeros(pair_size,1);
        eps_alpha = 0.01;
        for j = 1 : pair_size
            %         j = 154
            flag_agw = true(pair_size, 1);
            
            flag_agw(j) = false;
            tmp_index = index_pair(flag_agw);
            t1 = (pairs(j,5) - pairs(flag_agw,5)).^2 + (pairs(j,6) - pairs(flag_agw,6)).^2 + (pairs(j,7) - pairs(flag_agw,7)).^2;
            [value, agw_index] = sort(t1);
            u_vector = pairs(tmp_index(agw_index(1:neighbor_num)), 8:10);
            %         tmp_index(agw_index(1:neighbor_num))
            delta = 1.24 * mean(sqrt(value(1:neighbor_num)));
            weight = exp(-value(1:neighbor_num).^2/delta^2);
            sumn0 = sum(weight);
            
            agw_weight = weight./sum(weight);
            u_pos = agw_weight' * u_vector;
            
            u_median = median(u_vector);
            %         u_median = u_pos;
            rm = median(vecnorm(u_vector - u_median,2,2));
            u_cur = pairs(j,8:10);
            r0 = vecnorm(u_cur - u_median)/(rm + eps_alpha);
            
            %%%
            u_median_new = u_vector./(value(1:neighbor_num) + eps_alpha);
            r0_new = vecnorm( u_cur/(median(value(1:neighbor_num)) + eps_alpha) - median(u_median_new))/...
                (median(vecnorm(u_median_new - median(u_median_new),2,2)) + eps_alpha);
            rvalue(j) = r0;
        end
        index_col = index_agw(rvalue>tol);
        index_rank = assignment(index_agw(rvalue>tol));
        relapce_num = sum(rvalue > tol);
        for k = 1 : relapce_num
            %         dist_hugs(index_col(k), index_rank(k)) = inf;
            dist_hugs(index_rank(k), index_col(k)) = inf;
        end
        
        
    end
    parapairs{frame_idx} = pairs;
end

fprintf('\n');

%% boundary and mapping
res = x_res(2)-x_res(1);
boundaries = boundaries + 10;
boundary = [min(boundaries(:,1))-10, min(boundaries(:,2))-10,max(boundaries(:,3))-10, max(boundaries(:,4))-10];
new_res = 6*res;
x_res1 = boundary(1) : new_res : boundary(3);
y_res1 = boundary(2) : new_res : boundary(4);
z_res1 = z_res(1) : new_res : z_res(end);
size_x = length(x_res1);
size_y = length(y_res1);
size_z = length(z_res1);
size_x = floor(size_x/8)*8;
size_y = floor(size_y/8)*8;
size_z = floor(size_z/8)*8;
x_res2 = linspace(boundary(1),boundary(3),size_x);
y_res2 = linspace(boundary(2),boundary(4),size_y);
z_res2 = linspace(z_res(1),z_res(end),size_z);
new_res = min([x_res2(2)-x_res2(1),y_res2(2)-y_res2(1),z_res2(2)-z_res2(1)]);
x_res = x_res2;
y_res = boundary(2) : new_res : boundary(4);
z_res = z_res(1) : new_res : z_res(end);
size_x = floor(size_x/8)*8;
size_y = floor(size_y/8)*8;
size_z = floor(size_z/8)*8;
x_res = x_res(1:size_x);
y_res = y_res(1:size_y);
z_res = z_res(1:size_z);

x_res3 = 1 : length(x_res);
y_res3 = 1 : length(y_res);
z_res3 = 1 : length(z_res);
size_xy = size_x * size_y;
size_xyz = size_x * size_y * size_z;
npixels = size_xyz;
volume2d = zeros(size_xy,2);
for i = 1 : size_x
    volume2d((i-1)* size_y +1 : i * size_y,:) = [ones(size_y,1) * x_res(i),y_res'];
    volume2d_pos((i-1)* size_y +1 : i * size_y,:) = [ones(size_y,1) * x_res3(i),y_res3'];
end
volume3d = zeros(size_x * size_y * size_z, 3);
for i = 1 : size_z
    volume3d((i-1) * size_x * size_y+1 : i * size_x * size_y, :) = [volume2d, ones(size_x * size_y, 1) * z_res(i)];
    volume3d_pos((i-1) * size_x * size_y+1 : i * size_x * size_y, :) = [volume2d_pos, ones(size_x * size_y, 1) * z_res3(i)];
end

vol_left_current_map = zeros(size_xyz, 2, frame_num);
vol_right_current_map = zeros(size_xyz,2,frame_num);
occupancy_left_current = zeros(size_xyz,  frame_num);
pointnum=[];
for frame = frame_start : frame_stop
    frame_idx = frame - frame_start+1;
    pairs = parapairs{frame_idx};
    pointnum=[pointnum,size(pairs,1)];
    pairs1 = pairs(pairs(:,5)>=x_res(1) & pairs(:,5) < x_res(end) & ...
        pairs(:,6) >=y_res(1) & pairs(:,6) < y_res(end) & ...
        pairs(:,7) >=z_res(1) & pairs(:,7) < z_res(end),:);
    index=(round((pairs1(:,5)-x_res(1))/new_res))*size_y + round((pairs1(:,6)-y_res(1))/new_res) + 1 + round((pairs1(:,7)-z_res(1))/new_res) * size_xy;
    %     table1 = tabulate(index);
    index1 = 1 : size(pairs1,1);
    index2=index1(index<size_xyz);
    vol_left_current_map( index(index2), :, frame_idx) = pairs1(index2,12:13);
    vol_right_current_map(index(index2),:,frame_idx) = pairs1(index2,14:15);
    occupancy_left_current(index(index2),frame_idx) = 1;
end


pyramid_vol{1}.vol_left_current_map = vol_left_current_map;
pyramid_vol{1}.vol_right_current_map = vol_right_current_map;
pyramid_vol{1}.occupancy_left_current = occupancy_left_current;


