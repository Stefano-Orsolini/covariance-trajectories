% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Plot realtime covariance matrix as trajectories of
% coordinates (autocov_A, autocov_B, crosscov_AB).
% AB is any unique combination of sources from an 
% input signals matrix of size samples-by-sources.
%
% stefano.orsolini@gmail.com
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% analysis time window in seconds
params.t_show = 2.0;

% filter: 0 none, 1 global component removal, 2 bandpass filter
% NOTE: these settings are exclusive, cascading is not implemented
params.do_filt = 1;

% show crosscovariance boundary planes
params.show_plane = true;

% show indices of signals pairs
params.show_details = false;

% show signals timeseries (in another figure)
params.show_timeseries = false;

% ~~~ less relevant settings ~~~
% set plot limits
if params.do_filt==0
    params.ac_lim = 3500;
elseif params.do_filt==1
    params.ac_lim = 500;
elseif params.do_filt==2
    params.ac_lim = 1000;
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
    params.events =  importdata('S_events.dat');
end
% signals matrix: samples-by-sources
load('data/S_signals.mat');
params.signals = S_signals;
% sampling frequency
params.Fs = 160;

% IIR filter parameters
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
RTA_visualizer(params);
