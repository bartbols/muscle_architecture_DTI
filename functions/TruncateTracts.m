function varargout = TruncateTracts( DTItracts,SurfModel,varargin)
%%TRUNCATETRACTS Truncates the fiber tracts so that they terminate inside
% the muscle surface model and, optionally, outside the aponeurosis. The
% truncated length and fibre indices are returned.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% [ fibindex_trunc,length_trunc ] = TruncateTracts( DTItracts,surface)
% 
% ----------------- INPUT -----------------
% - DTItracts    :  a structure array containing at least the fields tracts_xyz
%                   and fibindex.
% - SurfModel    :  a structure array containing the fields 'vertices' and 
%                   'faces' of the surface model.
%
% Optional inputs, provided as 'parameter',<value> pairs:
% - aponeurosis   : a structure array containing the fields 'vertices' and 
%                   'faces' of the aponeurosis. A filename of an STL file
%                   can also be provided.
%         
% ----------------- OUTPUT -----------------
% - fibindex_trunc: n x 2 array with the first and last index of the
%    truncated fibre tracts that are within the surface
% - length_trunc: n x 1 array with the length of the truncated fibre tracts
%
% If one output argument is provided, the fields 'fibindex_trunc' and
% 'length_trunc' are added to DTItracts. If two outputs are provided, the
% first output is 'fibindex_trunc' and the second output 'length_trunc'
%
% Uses the function inpolyhedron.m to determine whether a point is inside
% or outside the surface.

%% Check inputs

tic
p = inputParser;
addRequired(p,'DTItracts',@(x) isstruct(x) || exist(x,'file')==2)
addRequired(p,'SurfModel',@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'aponeurosis',[])
parse(p,DTItracts,SurfModel,varargin{:})

aponeurosis = p.Results.aponeurosis;

%% Read inputs
% if tract filename is provided, read the file.
if ~isstruct(DTItracts)
    DTItracts = load(DTItracts);
end

% if surface model filename is provided, read the stl file
if ~isstruct(SurfModel)
    if exist(SurfModel,'file') == 2
        SurfModel = stlread(SurfModel);
    else
        error('%s does not exist.',SurfModel)
    end
end

% Read the aponeurosis surface, if provided
if ~isempty(aponeurosis) && ~isstruct(aponeurosis)
    if exist(aponeurosis,'file') == 2
        aponeurosis = stlread(aponeurosis);
    else
        error('%s does not exist.',aponeurosis)
    end
end

% Check if all required fields are available in DTItracts and surface
if ~isfield(DTItracts,'tracts_xyz')
    error('Required field tracts_xyz not found in DTItracts.')
end
if ~isfield(DTItracts,'fibindex')
    error('Required field fibindex not found in DTItracts.')
end

if ~isfield(SurfModel,'vertices')
    error('Required field vertices not found in SurfModel.')
end

if ~isfield(SurfModel,'faces')
    error('Required field faces not found in SurfModel.')
end

%% Fibre truncation
tic
nFib = size(DTItracts.fibindex,1);
fibindex_trunc = NaN(nFib,2);
length_trunc   = NaN(nFib,1);
hwait = waitbar(0,'Calculating which tract points are inside the muscle',...
    'Name','Progress bar TruncateTracts');

% Calculate points are inside the surface model. It is much faster to do 
% this for all points at once then inside the loop for each fibre 
% individually.
% INSIDE_muscle =  inside_surface(SurfModel,DTItracts.tracts_xyz');
INSIDE_muscle =  inpolyhedron(SurfModel,DTItracts.tracts_xyz');

if isempty(aponeurosis)
    % No aponeurosis is provide, so all points are outside the aponeurosis.
    OUTSIDE_apo = true(size(INSIDE_muscle));
else
    % An aponeurosis has been provided. Calculate which points are outside 
    % the aponeurosis.    
    waitbar(0,hwait,'Calculating which tract points are inside the aponeurosis')
%     OUTSIDE_apo = ~inside_surface(aponeurosis,DTItracts.tracts_xyz');
    OUTSIDE_apo = ~inpolyhedron(aponeurosis,DTItracts.tracts_xyz');
end
waitbar(0,hwait,sprintf('Truncating %d fibres',nFib))
for fibnr =  1:nFib
    if any(round(linspace(1,nFib,11))==fibnr)
        waitbar(fibnr/nFib,hwait)
    end
    
    % Remove tract points outside the muscle volume
    if DTItracts.fibindex(fibnr,1) < DTItracts.fibindex(fibnr,2)
        sgn = 1;
    else
        sgn = -1;
    end
    p = DTItracts.fibindex(fibnr,1) : sgn : DTItracts.fibindex(fibnr,2);
    
    % Check which points of this fibre are inside the muscle and
    % outside the aponeurosis
    dsig = diff([1 ~(INSIDE_muscle(p) & OUTSIDE_apo(p))' 1]);
    startIndex = find(dsig < 0);
    endIndex   = find(dsig > 0)-1;
    len = endIndex - startIndex+1;
    [nSteps,maxIndex] = max(len);

%   % Plot the selected points along the tract
%     pts = DTItracts.tracts_xyz(:,p);
%     plot3(pts(1,startIndex:endIndex),...
%           pts(2,startIndex:endIndex),...
%           pts(3,startIndex:endIndex),...
%         'go','MarkerSize',6)
%     
%     plot3(pts(1,:),...
%           pts(2,:),...
%           pts(3,:),...
%         'yo','MarkerSize',6)
    
    if isempty(len) || nSteps < 5
        % No or too few points are inside the surface and outside the
        % aponeurosis
        continue
    else
        % Select the longest continuous section inside the muscle surface
        fibindex_trunc(fibnr,1) = p(startIndex(maxIndex));
        fibindex_trunc(fibnr,2) = p(endIndex(maxIndex));

        % Calculate stepsize as the distance between the first and
        % second tract points. It is assumed to fibre tracking was
        % performed with a constant stepsize.
        stepsize = norm(DTItracts.tracts_xyz(1:3,2)-DTItracts.tracts_xyz(1:3,1));
        length_trunc(fibnr)     = (nSteps-1) * stepsize;
    end
    % Make sure the startpoint (column 1 in fibindex_trunc) has a lower z-value than
    % the endpoint (column2 in fibindex_trunc).
    if DTItracts.tracts_xyz(3,fibindex_trunc(fibnr,1)) > DTItracts.tracts_xyz(3,fibindex_trunc(fibnr,2))
        fibindex_trunc(fibnr,:) = fliplr(fibindex_trunc(fibnr,:));
        DTItracts.fibindex(fibnr,:) = fliplr(DTItracts.fibindex(fibnr,:));
        
    end
end

% Return as separate output arguments if 2 outputs are requested or add to
% structure 'DTItracts' when one output is requested.
if nargout == 1
    DTItracts.fibindex_trunc = fibindex_trunc;
    DTItracts.length_trunc   = length_trunc;
    varargout{1} = DTItracts;
elseif nargout == 2
    varargout{1} = fibindex_trunc;
    varargout{2} = length_trunc;
end
t_elapsed = toc;
fprintf('It took %.2f seconds to truncate fibre tracts.\n',t_elapsed)
% fprintf('completed\n')
close(hwait)
end % of function

