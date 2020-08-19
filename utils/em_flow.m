function [flow_out, time_shifted_points_out, event_marker_out] = em_flow(...
    events, ...
    feature_pos, ...
    flow_init, ...
    params)
%% em1 flow
flow = flow_init;
flow_out = [nan; nan];
time_shifted_points_out = [];
delta_flows = zeros(params.em1_params.max_iters, 1);
num_iter = 0;
event_window = [];
n_in_window = 0;
target_time = events(1,1);
event_marker_out = [];
%%% translate the events so that they are centered at the last
%%% position of the feature.
centered_events = events;
centered_events(:, 2:3) = centered_events(:, 2:3) - feature_pos;
prev_flow = flow;

while true
    if num_iter > params.em1_params.max_iters
        return;
    end
    %%% shift the centered events in time by the current
    %%% estimate of the flow
    time_shifted_points = centered_events(:, 2:3) + flow' .* (target_time - centered_events(:,1));
    
    %%% event window computing
    if isempty(event_window)
        event_window = time_shifted_points(:, 1) >= - params.window_size/2 & ...
            time_shifted_points(:,2) >= -params.window_size/2 & ...
            time_shifted_points(:,1) <= params.window_size/2 & ...
            time_shifted_points(:,2) <= params.window_size/2;
        n_in_window = sum(event_window);
        %%% ignore if there are not enough events
        if n_in_window < params.min_events_for_em
            return
        end
        centered_events1 = centered_events(event_window, :);
        time_shifted_points1 = time_shifted_points(event_window, :);
    end
    %%% scale the events so that the computed distances are
    %%% equal to (x1 - x2)/(2 * sigma^2). Saves us having to
    %%% divide all the correspondences later.
    normalized_events = time_shifted_points1 /(sqrt(2) * params.em1_params.sigma);
    
    %%% The KD tree allows us to compute distances between
    %%% neighboring events, while also performing outlier
    %%% rejection. Unfortunately, as the time shifted events
    %%% change every iteration, it must also be reconstructed
    %%% at every iteration.
    
    kdtree = KDTreeSearcher(...
        normalized_events, ...
        'Distance',...
        'euclidean');
    %%% Threshold distance by the Malhalanobis distance
    [neighbors_cell, distance_cell] = rangesearch(...,
        kdtree, ...
        normalized_events, ...
        params.em1_params.max_distance);
    
    distancestacked = cell2mat(cellfun(...
        @transpose, distance_cell, 'UniformOutput', false))';
    
    %%% Can't solve for flow without any correspondences
    if isempty(distancestacked)
        return
    end
    %%% Number of neighbors for each event
    %%% Note: 'length' is not the same as @length
    num_neighbors_per_event = cellfun('length', neighbors_cell);
    neighbor_inds = cell2mat(cellfun(@transpose, neighbors_cell, 'UniformOutput', false))';
    event_inds = repelem(1 : n_in_window, num_neighbors_per_event);
    
    %%% The event to neighbor graph is undirected, so no need
    %%% to double count the event-neighbor correspondences.
    
    valid_correspondences = neighbor_inds > event_inds;
    neighbor_inds = neighbor_inds(valid_correspondences);
    event_inds = event_inds(valid_correspondences);
    distancestacked = distancestacked(valid_correspondences);
    
    %%% Distance are already scaled by the variance. Note that
    %%% the original equation in the paper is the sum of two
    %%% weights. Here we simplify it with a single weight for
    %%% speed.
    
    weights = exp(-distancestacked);
    %%% Simultaneously minimize over flow and translation
    neighbor_events = centered_events1(neighbor_inds, : );
    original_events = centered_events1(event_inds, :);
    
    %%% These are the X and D matrices specified in the
    %%% original paper, except the weighting is handled by
    %%% multiplying the D with the full weight.
    X = original_events(: , 2 : 3) - neighbor_events(: , 2 : 3);
    D = original_events(:, 1) - neighbor_events(:, 1);
    weighted_D = bsxfun(@times, D , weights');
    DDT = weighted_D' * D ;
    XDT = weighted_D' * X;
    %                 DDT = D;
    %                 XDT = X;
    
    flow = XDT'/DDT';
    
    %%% Calculate change in flow, plot debug information
    if( norm(flow - prev_flow) < params.em1_params.min_err)
        break;
    end
    delta_flow(num_iter + 1 ) = norm(flow - prev_flow);
    prev_flow = flow;
    num_iter = num_iter + 1;
end


%%% Calculate the final shifted events to be used in EM2 later.
dt = events(end, 1) - events(1, 1);
centered_events = events;
centered_events(:, 2:3) = centered_events(:, 2:3) - feature_pos + flow' * dt;
time_shifted_points2 = centered_events(:, 2: 3) + bsxfun(@times, flow', (events(end, 1) - centered_events(:, 1)));
%%% Make the window a little bigger.
window_size = round(params.window_size * 1.5);
event_window = time_shifted_points2(:, 1) >= -window_size/2 & ...
    time_shifted_points2(:, 2) >= - window_size/2 & ...
    time_shifted_points2(:, 1) <= window_size/2 & ...
    time_shifted_points2(:, 2) <= window_size/2;
time_shifted_points_out = time_shifted_points2(event_window, :);

%%% Assign the event marker
window_size2 = round(params.window_size * params.window_mag);
event_window = time_shifted_points2(:, 1) >= -window_size2/2 & ...
    time_shifted_points2(:, 2) >= - window_size2/2 & ...
    time_shifted_points2(:, 1) <= window_size2/2 & ...
    time_shifted_points2(:, 2) <= window_size2/2;
event_marker_out = event_window;
flow_out = flow;

end