function [ DTItracts ] = CalcArchitecture( DTItracts,SurfModel,varargin)
%CALCARCHITECTURE truncates fibre tracts to the muscle surface model, then
% fits polynomials to the tracts and extrapolates them to their
% intersection with the muscle surface. Fibre lengths, pennation angles and
% curvatures are calculated for all tracts.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% August 2018
%
% Please refer to the following paper when using this script:
%
% Bolsterlee B, D'Souza A, Gandevia SC, Herbert RD.
% How does passive lengthening change the architecture of the human medial gastrocnemius muscle?
% J Appl Physiol. 2017;122(4):727-38.
%
% ----------------- USAGE ----------------- 
% DTItracts = CalcArchitecture(DTItracts,SurfModel)
% or
% DTItracts = CalcArchitecture(DTItracts,SurfModel,order)
%
% ----------------- INPUT -----------------
% Required:
% - DTItracts : a structure array containing at least the fields 'tracts_xyz' 
%               and 'fibindex' OR the filename of a MAT-file containing
%               these variables.
% 
% - SurfModel : a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the muscle 
%               surface used for extrapolating the tracts OR the filename
%               of an STL-file containing this surface.
%
% Optional inputs, provided as 'parameter',<value> pairs:
% - order     : order of the polynomial fit. Default = 3
% - aponeurosis: a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               aponeurosis used for extrapolating the tracts.
%
% tracts_xyz is a 3 x N array containing the xyz-coordinates of all tracts
% points of all fibres in the set.
% fibindex is an M x 2 array containing the start index (first column) and 
% end index (second column) of the tracts points of a fibre. The number of
% rows M equals the number of fibres in the set.
%
% Example:
%
% fibindex(1,:) = [2 120] means that fibre #1 has tracts_xyz(:,2) as
% start point, tracts_xyz(:,3) as second point, tracts_xyz(:,4) as third
% point, ... and tracts_xyz(:,120) as end point.
%
% fibindex(2,:) = [121 141] means that fibre #2 has tracts_xyz(:,121) as
% start point, tracts_xyz(:,122) as second point, tracts_xyz(:,123) as third
% point, ... and tracts_xyz(:,141) as end point.
% 
%
% ----------------- OUTPUT -----------------
% - DTItracts : a structure array with reconstruction parameters and the
%               architectural measures added as fields.
%
% This function calls TruncateTracts, ExtrapolateTracts, CalcPenAngle and
% CalcCurvature. Type help <function_name> for more information on the
% outputs of each of these functions.

% Optional: 
% - order     : order of the polynomial fit. Default = 3
%% Check inputs
p = inputParser;
addRequired(p,'DTItracts',@(x) isstruct(x) || exist(x,'file')==2)
addRequired(p,'SurfModel',@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'order',3,@(x) isscalar(x) && x>0)
addParameter(p,'aponeurosis',[])
parse(p,DTItracts,SurfModel,varargin{:})

order       = p.Results.order;
aponeurosis = p.Results.aponeurosis;

%% Read and check inputs
% If a filename is provided, read the file with the DTItracts.
if ~isstruct(DTItracts)
    if exist(DTItracts,'file') == 2
        DTItracts = load(DTItracts);
    else
        error('%s does not exist.',DTItracts)
    end
end

% Check if DTItracts contains variables tracts_xyz and fibindex.
if ~isfield(DTItracts,'tracts_xyz')
    error('DTItracts does not contain a field tracts_xyz, a 3xn array of tract point coordinates.')
end
if ~isfield(DTItracts,'fibindex')
    error('DTItracts does not contain a field fibindex, a m x 2 array with start and end indices of for each fibre.')
end


% If a filename is provided, read the stl file with the surface model.
if ~isstruct(SurfModel)
    if exist(SurfModel,'file') == 2
        SurfModel = stlread(SurfModel);
    else
        error('%s does not exist.',SurfModel)
    end
end

% Read the aponeurosis surface, if an aponeurosis is defined.
if ~isempty(aponeurosis) && ~isstruct(aponeurosis)
    if exist(aponeurosis,'file') == 2
        aponeurosis = stlread(aponeurosis);
    else
        error('%s does not exist.',aponeurosis)
    end
end
%%
% Truncate tracts
DTItracts   = TruncateTracts(DTItracts,SurfModel,...
    'aponeurosis',aponeurosis);

% Extrapolate tracts
DTItracts   = ExtrapolateTracts(DTItracts,SurfModel,...
    'order',order,...
    'aponeurosis',aponeurosis);

% Calculate pennation angle
DTItracts.penangle = CalcPenAngle( DTItracts,SurfModel,...
    'aponeurosis',aponeurosis);
    
% Calculate curvature
DTItracts.curvature = CalcCurvature( DTItracts.PolyCoeff);

% Add angle between direction vectors at endpoint. This angle can later be
% used to exclude fibres that originate and insert on the same aponeurosis.
DTItracts.ang = acosd(sum(squeeze(DTItracts.endpoints_dir(:,1,:)) .* squeeze(DTItracts.endpoints_dir(:,2,:)),2)); 

end

