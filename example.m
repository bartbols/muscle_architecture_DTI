clear
addpath(genpath('functions'))

%% Fascicle reconstruction WITHOUT an internal aponeurosis
% Set the filename of the fibre tracts and the muscle surface model.

% Fibre tracts have to be provided as a MAT-file containing at least the
% following variables:
% - tracts_xyz
% - fibindex
% Type 'help CalcArchitecture' for more information on the expected format 
% of these variables.
T = load('tracts_MG.mat');

% Set the filename of the surface model. If a filename is provided, the
% surface is expected in STL-format. Alternatively, a structure
% containing the variables 'vertices' and 'faces' (in Matlab's Patch Object
% format) can be provided.
M = 'MG.stl';

% Reconstruct the architecture by truncating fibre tracts to the boundary of
% the 3D surface model, fitting polynomials on the raw fibre tracts, and
% then extrapolating the polynomials to their intersection with the surface
% model. Architectural parameters are added as fields to the structure 'T'.

% Type 'help CalcArchitecture in the command window for more information
% about optional inputs to CalcArchitecture.
T = CalcArchitecture(T,M);

% Save the extended fibre tracts and muscle architecture measurements as a 
% new file.
save('tracts_MG_ext.mat','-struct','T')

% Inspect the results with the FibreInspecter. Type 'FibreInspecter' in the
% command prompt and load the tract and surface files':
% To load the tracts, click on Load files > Load tracts... in the menu bar
% and select 'tracts_MG_ext.mat'.
% To load the surface model, click on Load files > Load surface... in the 
% menu bar and and select 'MG.stl')

%% Fascicle reconstruction WITH an internal aponeurosis
% Set the filename of the fibre tract data, the muscle surface model and
% the aponeurosis surface model.
clear
T = load('tracts_TA.mat');
M = 'TA.stl';
A = 'TA_apo.stl';

% Reconstruct the architecture by truncating fascicles to the boundary of
% the 3D surface model or the aponeurosis, fitting polynomials on the raw
% fibre tracts, and then extrapolating the polynomials to their
% intersection with the surface/aponeurosis model. Architectural parameters 
% are added as fields to the structure 'T'.
T = CalcArchitecture(T,M,'aponeurosis',A);

% Save the extended fibre tracts and muscle architecture measurements as a 
% new file.
save('tracts_TA_ext.mat','-struct','T')

% Inspect the results again with the FibreInspecter. Now also load the
% aponeurosis model.
