% clear all;
% clc;
%% PIV 2d reconstruction with high speed
% maximum iterations of EM
em1_params.max_iters = 50;
% Convariance of traget point Gaussians
em1_params.sigma = 2;
% min change in error before optimization stops
em1_params.min_err = 1.0;
% Max distance for rangesearch for em1 beyond which points are considered
% outliers
em1_params.max_distance = 1.0;
em1_params.centered = 4;
params.em1_params = em1_params;
% window size
params.window_size = 15;
% window search size
params.window_mag = 1.2;
% minimum events for em
params.min_events_for_em = 20;
% max events in one time window
params.max_events_per_window = 100000;
params.min_distance = 4;