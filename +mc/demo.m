%% Setting up experiment
clear;clc;
target.Dims     = 2;                        % number of sampling dimensions
target.mu       = 5*ones(1,target.Dims);    % initial starting point or mean location
target.Sigma    = eye(target.Dims);         % initial covariance
[target.Sigma,target.chol,~] = make_symm(target.Sigma);
target.M        = 1e3*target.Dims;          % number of parallel samples/chains
target.name     = 'bimodal';                % name of target
target.log      = @mc.target.bimodal;       % function handle for target log-density 
target.inst     = @(x) log_mvnpdf(x,target.mu,target.Sigma);    % function handle for instrumental log-density 
target.rand     = @(target) q_mvnrnd(target);   % random number generator for instrumental density

%% Example: Importance sampling
    [particles] = mc.is.sampling( target );
%% Example: Adaptive importance sampling
    [target,particles] = mc.is.adaptive_sampling( target );
%% Example: Sequential Monte Carlo
particles = mc.smc.sampling( target );
%% Example: Slice sampling 
% Working out the number of files to store the samples, each file is at
% lease approximately 100MB
target.m    = ceil(1e8/filesize(particles(end))); 
target.rep  = 100; % Desired number of slice sampling repetitions
target.div  = ceil(target.rep/target.m); % number of files
particles   = mc.mc.slice.sampling( target, particles );
figure(1);hist3(mc.reproc.aggregate(particles),50*[1 1]);
%% Adaptive Metropolis
target.dsp = 1;     % Display status of burnin runs, 1 - display | 0 - suppress
target.rep = 1e5;   % Desired number of Metropolis sampling repetitions
target.mod = 0;     % Mode for generating auxiliary innovation
target     = mc.mc.rwm.adapt(target);
particles  = mc.mc.rwm.sampling(target);
figure(2);hist3(mc.reproc.aggregate(particles),50*[1 1]);

%% Cleaning up folder 
delete *.mat; clear; close all;
