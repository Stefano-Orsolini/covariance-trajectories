% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot realtime covariance matrix as trajectories of
% coordinates (autocov_A, autocov_B, crosscov_AB).
% AB is any unique combination of sources from an 
% input signals matrix of size samples-by-sources.
%
% dependancies:
% - FilterM (File Exchange)
% - designfilt (Signal Processing Toolbox)
%
% stefano.orsolini@gmail.com
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% analysis time window in seconds
params.t_show = 2.0;
% apply global component removal
params.do_global = true;
% apply bandpass filter
params.do_filter = false;
% show crosscovariance boundary planes
params.show_plane = true;
% show indices of signals pairs
params.show_details = false;
% show signals timeseries (in another figure)
params.show_timeseries = false;

% ~~~ less relevant settings ~~~
% set plot limits
if params.do_global
    if params.do_filter
        params.ac_lim = 100;
    else
        params.ac_lim = 500;
    end
else
    if params.do_filter
        params.ac_lim = 1000;
    else
        params.ac_lim = 3500;
    end
end


% fix axes on limits
params.fixed_axes = true;

% zoom locally on last samples
params.show_local = false;
% frames per second
params.FPS = 25;
% visualization speed multiplier
params.timescale = 1;
% show experimental paradigm label textbox
params.show_box = true;
% loop visualization
params.do_loop = true;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if params.show_box
    % experimental paradigm sequence
    params.events =  importdata(fullfile('data','S_events.dat'));
end
% signals matrix: samples-by-sources
load(fullfile('data','S_signals.mat'));
params.signals = S_signals;
% sampling frequency
params.Fs = 160;

% bandpass filter parameters
d1 = designfilt('bandpassiir', ...
    'StopbandFrequency1', 5, ...
    'PassbandFrequency1', 7, ...
    'PassbandFrequency2', 35, ...
    'StopbandFrequency2', 37, ...
    'StopbandAttenuation1', 6, ...
    'PassbandRipple', 1, ...
    'StopbandAttenuation2', 6, ...
    'SampleRate', params.Fs, ...
    'MatchExactly', 'passband');

% % visual check the filter
%fvtool(d1)

% pass transfer function coefficients
[params.b, params.a] = tf(d1);

% function call
viewer_covtraj(params);
