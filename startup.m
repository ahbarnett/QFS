function startup(dim)
% QFS matlab initialization.
% usage: startup or startup(2) does 2D QFS; startup(3) does 3D.
% Crucial to get your ambient dimension right then stick to it!
if nargin<1, dim=2; end

% path
h = fileparts(mfilename('fullpath'));        % direc of this file
addpath(h);
addpath([h '/utils']);

% other utils needed
addpath ~/numerics/linalg/randutv/matlab

% my matlab prefs
set(0,'showHiddenHandles','on');        % eg lines linesmoothing
format long g
format compact
set(groot,'DefaultAxesColorOrder',[0 0 1; 0 .5 0; 1 0 0; 0 .75 .75; ...
                    .75 0 .75; .75 .75 0; .25 .25 .25])
set(groot,'DefaultFigureColormap',jet(256))
set(groot,'defaultLegendAutoUpdate','off')

bie2d = '~/BIE2D';         % user to edit based on their installation dirs
bie3d = '~/BIE3D';

% cleanly access only one package...
warning('off','MATLAB:rmpath:DirNotFound');
if dim==2
  rmpath(genpath(bie3d));
  addpath(bie2d);  
  bie2dsetup
  addpath([h '/2D']);
elseif dim==3
  rmpath(genpath(bie2d));
  addpath(bie3d);
  bie3dsetup
  rmpath(genpath([bie3d '/timedomainwaveeqn']));    % gateway clash w/ randutv
  addpath([h '/3D']);
end
