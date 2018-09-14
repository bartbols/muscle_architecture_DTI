function penangle = CalcPenAngle( DTItracts,surf_model,varargin )
%CALCPENANGLE3D calculates the pennation angle of fibres described by the
%polynomial coefficient in PolyCoeff and endpoints, which are both fields in
%DTItracts, with the surface model surf_model, described by its fields
%'vertices' and 'faces'
% Pennation angle is defined as 90 degr minus the median angle between the 
% fibre and the normal vectors of all faces of the muscle surface/
% aponeurosis model within the search radius around the endpoint (default 
% value is 1.5mm).
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE -----------------
% DTItracts = TruncateTracts( DTItracts,surf_model)
% or
% DTItracts = TruncateTracts( DTItracts,surf_model,radius)
%
% ----------------- INPUT -----------------
% Required
% - DTItracts  : a structure array containing at least the fields PolyCoeff
%               and endpoints, as generated with ExtrapolateTracts
% - surf_model : a structure array 'surface' containing the fields
%                'vertices' and 'faces' and optionally 'normals' or a
%                filename of an STL-file of the surface model.
%
% Optional input, provided as 'parameter',<value> pairs:
% - aponeurosis   : a structure array containing the fields 'vertices' and
%                   'faces' of the aponeurosis. A filename of an STL file
%                   can also be provided.
% - radius    : the radius around the fibre endpoint in which triangles are
%               included for pennation angle calculation. Default = 2.5
%
% ----------------- OUTPUT -----------------
% penangle    : n x 2 array of pennation angle of endpoint 1 and 2. n is
%               the number of fibres

%% Check inputs
p = inputParser;
addRequired(p,'DTItracts')
addRequired(p,'surf_model')
addParameter(p,'radius',1.5,@(x) validateattributes(x,{'numeric'},{'scalar'}) )
addParameter(p,'aponeurosis',[])

parse(p,DTItracts,surf_model,varargin{:})
radius      = p.Results.radius;
aponeurosis = p.Results.aponeurosis;

% Read the surface model if a filename is provided.
if ~isstruct(surf_model)
    surf_model = stlread(surf_model);
end

% Read the aponeurosis surface, if provided
if ~isempty(aponeurosis) && ~isstruct(aponeurosis)
    if exist(aponeurosis,'file') == 2
        aponeurosis = stlread(aponeurosis);
    else
        error('%s does not exist.',aponeurosis)
    end
end

% if tract filename is provided, read the file.
if ~isstruct(DTItracts)
    DTItracts = load(DTItracts);
end

% Check if pennation angle is defined relative to the aponeurosis. If so,
% check if the aponeurosis is provided.
if any(DTItracts.attach_type(:)==2) && isempty(aponeurosis)
    error('Pennation angle is defined relative to the aponeurosis, but the aponeurosis model is not provided.')
end

%%
% Calculate the centre of all faces of the muscle model
tic
C_mus = squeeze(mean(reshape(surf_model.vertices(surf_model.faces,:)',3,size(surf_model.faces,1),3),3))';

if ~isempty(aponeurosis)
    % Calculate the centre of all faces of the aponeurosis model
    C_apo = squeeze(mean(reshape(aponeurosis.vertices(aponeurosis.faces,:)',3,size(aponeurosis.faces,1),3),3))';
    % Calculate face normals of aponeurosis model if not provided
    aponeurosis.normals = facenormals(aponeurosis);
    
end

% Calculate face normals of muscle surface model
surf_model.normals = facenormals(surf_model);

nFibres = size(DTItracts.endpoints,1);
penangle = NaN(nFibres,2);
for fibnr = 1:nFibres
    for ep = 1:2
        
        % Select triangles within the selected radius around the end point
        if DTItracts.attach_type(fibnr,ep) == 1
            % ... if attachment is on the muscle surface
            C = C_mus;
            model = surf_model;
        elseif DTItracts.attach_type(fibnr,ep) == 2
            % ... if attachment is on the aponeurosis
            C = C_apo;
            model = aponeurosis;
        else
            continue
        end
        tmp = C - repmat(squeeze(DTItracts.endpoints(fibnr,ep,:))',size(C,1),1);
        dist2centre = sqrt(sum(tmp.^2,2));
        idx = find(dist2centre < radius);
        
        %         if ep == 1
        %             t = DTItracts.PolyCoeff(fibnr).t0;
        %         else
        %             t = DTItracts.PolyCoeff(fibnr).t1;
        %         end
        %         slope = [polyval(polyder(DTItracts.PolyCoeff(fibnr).x),t) ,...
        %             polyval(polyder(DTItracts.PolyCoeff(fibnr).y),t) ,...
        %             polyval(polyder(DTItracts.PolyCoeff(fibnr).z),t)];
        %         slope = slope/norm(slope); % make unit vector
        
        slope = squeeze(DTItracts.endpoints_dir(fibnr,ep,:))';
        
        
        % Only include vectors that are pointing in the opposite direction
        % to the slope at the fascicle endpoint.
        proj = sum(model.normals(idx,:) .* repmat(slope,length(idx),1),2);
        if DTItracts.attach_type(fibnr,ep) == 1
            idx(proj<0) = [];
        elseif DTItracts.attach_type(fibnr,ep) == 2
            idx(proj>0) = [];
        end
        
        % calculate angle with normal vectors of selected triangles
        if isempty(idx)
            continue
        end
%         figure
%         InspectTracts('Tracts',DTItracts,...
%             'SurfModel',surf_model,...
%             'aponeurosis',aponeurosis,...
%             'Selection',fibnr,...
%             'ToPlot','poly',...
%             'PlotStats',false,...
%             'LineWidth',2)
%         hold on
%         
%         patch('Vertices',model.vertices,...
%             'Faces',model.faces(idx,:),...
%             'FaceColor','g',...
%             'FaceAlpha',0.2,...
%             'EdgeColor','g',...
%             'LineWidth',2)
%         quiver3(C(idx,1),C(idx,2),C(idx,3),...
%             model.normals(idx,1),...
%             model.normals(idx,2),...
%             model.normals(idx,3),5,'r')
%         
%         quiver3(DTItracts.endpoints(fibnr,ep,1),...
%             DTItracts.endpoints(fibnr,ep,2),...
%             DTItracts.endpoints(fibnr,ep,3),...
%             slope(1),...
%             slope(2),...
%             slope(3),5,'k','LineWidth',3)
%         light
%         axis equal
%         
        
        penangle(fibnr,ep) = abs(median(asind(sum(model.normals(idx,:) .* repmat(slope,length(idx),1),2))));
        
        
    end
end
t_elapsed = toc;
fprintf('It took %.2f seconds to calculate pennation angles.\n',t_elapsed)
end % of function
