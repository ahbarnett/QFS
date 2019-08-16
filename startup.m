% QFS matlab initialization

% path
h = fileparts(mfilename('fullpath'));        % direc of this file
addpath(h);
addpath([h '/utils']);

% my matlab prefs
set(0,'showHiddenHandles','on');        % eg lines linesmoothing
format long g
format compact
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                    .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
set(groot,'defaultLegendAutoUpdate','off')

% 2D
addpath ~/BIE2D
bie2dsetup
