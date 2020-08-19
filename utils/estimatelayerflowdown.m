function uv = estimatelayerflowdown(position, uv, lambda1,  lambda2 , lambda3 , ...
    INNER_ITER, maxwarping, x2, y2, z2, vol_master_current_map, vol_slave_current_map , trans1, trans2, frame_num, levels,pyramid_levels)
[H, w, D] = size(position(:,:,:,1));
% frame_num = size(position,4);
npixels = H * w *D;

dx_kernel = zeros(3,3,3);
dy_kernel = zeros(3,3,3);
dz_kernel = zeros(3,3,3);
dx_kernel(2,2,2) = -1;
dx_kernel(3,2,2) =  1;
dy_kernel(2,2,2) = -1;
dy_kernel(2,3,2) = 1;
dz_kernel(2,2,2) = -1;
dz_kernel(2,2,3) = 1;

lap_x = zeros(3,3,3);
lap_x(2,2,1) = -1;
lap_x(2,2,2) = 6;
lap_x(1,2,2) = -1;
lap_x(2,1,2) = -1;
lap_x(2,3,2) = -1;
lap_x(3,2,2) = -1;
lap_x(2,2,3) = -1;
param_test = 0.4;
lambda1 = lambda1*param_test;
lambda2 = lambda2*param_test;
lambda3 = lambda3*param_test;
bound_cond = 'replicate'; %% or symmetric


vol_master_current_map = reshape(vol_master_current_map,[npixels,2,frame_num]);
vol_slave_current_map = reshape(vol_slave_current_map,[npixels,2,frame_num]);
for iteout = 1 : INNER_ITER
    fprintf('outer iteration %d \n', iteout);
    fprintf('frame processed: ');
    for i = 1 : frame_num
        fprintf('%d..',i);
        v1tx = vol_master_current_map(:, 1, i);
        v1ty = vol_master_current_map(:, 2, i);
        v2tx = vol_slave_current_map(:, 1,i);
        v2ty = vol_slave_current_map(:, 2,i);
            for ite = 1 : maxwarping
                %% This is to solve the Ax = b problem
                [A,b] = constructAmatrixdown(position, dx_kernel, dy_kernel, dz_kernel, lap_x, uv, x2, y2, z2, bound_cond, ...
                    param_test, lambda1, lambda2 , lambda3, npixels, iteout, i, frame_num, v1tx, v1ty, v2tx, v2ty, trans1, trans2, levels,pyramid_levels);
                
                [deltauv,~, ~, ~] = pcg(A,b,1e-3,100);
                deltauv(deltauv>1) = 1;
                deltauv(deltauv<-1) = -1;
                deltauv = reshape(deltauv, size(uv(:,:,:,:,1)));
                uv(:,:,:,:,i) = uv(:,:,:,:,i) + deltauv;
            end
    end
    fprintf('\n');
end
end
