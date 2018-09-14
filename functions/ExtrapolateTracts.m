function varargout = ExtrapolateTracts( DTItracts,SurfModel,varargin )
%EXTRAPOLATETRACTS This function fits a polynomial on the truncated fibre
% tracts and then extrapolates this polynomial linearly at its endpoints
% until the surface is intersected.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% ----------------- USAGE ----------------- 
% DTItracts = ExtrapolateTracts(DTItracts,SurfModel,order)
% or 
% [PolyCoeff,fibrelength,ext,pct_ext,endpoints,residual] = ExtrapolateTracts(DTItracts,SurfModel,order)
%
% ----------------- INPUT -----------------
% Required:
% - DTItracts : a structure array containing at least the fields tracts_xyz 
%               and fibindex_trunc (created with TruncateTracts.m).
% - SurfModel: a structure array containing the fields
%               'vertices','faces' with the vertices and faces of the 
%               surface used for extrapolating the tracts.
%
% Optional: 
% - order     : order of the polynomial fit. Default = 3
%
% ----------------- OUTPUT -----------------
% If one output argument is provided DTItracts is returned with the 
% following fields appended (n = number of fibres)
%
% - PolyCoeff   : n x 1 structure array with fields x, y and z containing the
%                 polynomial coefficients of the polynomial curve for the x,
%                 y and z coordinates, respectively. Also has fields t0
%                 and t1 with the first and last t-value of the curve.
% - fibrelength : n x 1 array with fibre length (curved length + extrapolated parts)
% - ext         : n x 2 array with extension (in mm) at each of the endpoints
% - pct_ext     : n x 1 array with total extension as percentage of the total fibre length
% - endpoints   : n x 2 x 3 array with xyz-location of extrapolated endpoints on the surface.
% - endpoints_dir: n x 2 x 3 array with xyz-vector of the slope at the
%                  endpoint
% - residual    : residual of polynomial fit (absolute mean distance of tracts
%                points to nearest point on the polynomial curve).
%
% If multiple output arguments are provided, the outputs are (in this
% order):
%
% [PolyCoeff,fibrelength,ext,pct_ext,endpoints,residual]
%

%% Check inputs
tic
p = inputParser;
addRequired(p,'DTItracts',@(x) isstruct(x) || exist(x,'file')==2)
addRequired(p,'SurfModel',@(x) isstruct(x) || endsWith(x,'.stl','IgnoreCase',true))
addParameter(p,'aponeurosis',[])
addParameter(p,'order',3,@(x) isscalar(x) && x>0)
parse(p,DTItracts,SurfModel,varargin{:})

order       = p.Results.order;
aponeurosis = p.Results.aponeurosis;

persistent ShowWarningMex
if isempty(ShowWarningMex);ShowWarningMex = true;end
if exist('FindNearestT_mex','file') == 3
    use_mex = true;
else
    use_mex = false;
    if ShowWarningMex == true
        warning('The MEX-file FindNearestT_mex is not available on the MATLAB path. Extrapolation may be slow.')
        ShowWarningMex = false; % only show the warning message once
    end
end

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

% Calculate edge vectors etc of muscle model once to speed up performance 
% of inside_surface and intersectRaySurface
% Get all vertices locations
data_musc.V0 = SurfModel.vertices(SurfModel.faces(:,1),:);
data_musc.V1 = SurfModel.vertices(SurfModel.faces(:,2),:);
data_musc.V2 = SurfModel.vertices(SurfModel.faces(:,3),:);

% Calculate edge vectors
data_musc.U = data_musc.V1 - data_musc.V0; 
data_musc.V = data_musc.V2 - data_musc.V0;

% Calculate dot product of edge vector
data_musc.UU = sum(data_musc.U.*data_musc.U,2);
data_musc.UV = sum(data_musc.U.*data_musc.V,2);
data_musc.VV = sum(data_musc.V.*data_musc.V,2);
data_musc.DD = data_musc.UV.*data_musc.UV - data_musc.UU.*data_musc.VV;

% Calculate normal of faces
data_musc.N = cross(data_musc.U,data_musc.V);


% Read the aponeurosis surface, if provided
if ~isempty(aponeurosis) && ~isstruct(aponeurosis)
    if exist(aponeurosis,'file') == 2
        aponeurosis = stlread(aponeurosis);
    else
        error('%s does not exist.',aponeurosis)
    end
end

if ~isempty(aponeurosis)
    % Calculate edge vectors etc of the aponeurosis model to speed up 
    % performance of inside_surface and intersectRaySurface
    
    % Get all vertices locations
    data_apo.V0 = aponeurosis.vertices(aponeurosis.faces(:,1),:);
    data_apo.V1 = aponeurosis.vertices(aponeurosis.faces(:,2),:);
    data_apo.V2 = aponeurosis.vertices(aponeurosis.faces(:,3),:);

    % Calculate edge vectors
    data_apo.U = data_apo.V1 - data_apo.V0; 
    data_apo.V = data_apo.V2 - data_apo.V0;

    % Calculate dot product of edge vector
    data_apo.UU = sum(data_apo.U.*data_apo.U,2);
    data_apo.UV = sum(data_apo.U.*data_apo.V,2);
    data_apo.VV = sum(data_apo.V.*data_apo.V,2);
    data_apo.DD = data_apo.UV.*data_apo.UV - data_apo.UU.*data_apo.VV;

    % Calculate normal of faces
    data_apo.N = cross(data_apo.U,data_apo.V);
end

% Check if all required fields are available in DTItracts and surface
if ~isfield(DTItracts,'tracts_xyz')
    error('Required field tracts_xyz not found in DTItracts.')
end
if ~isfield(DTItracts,'fibindex_trunc')
    error('Required field fibindex_trunc not found in DTItracts.')
end

if ~isfield(SurfModel,'vertices')
    error('Required field vertices not found in SurfModel.')
end

if ~isfield(SurfModel,'faces')
    error('Required field faces not found in SurfModel.')
end

%%
% Number of fibres
nFib = size(DTItracts.fibindex_trunc,1);

% Use the truncated tracts for fitting a polynomial and extrapolation.
if ~isfield(DTItracts,'fibindex_trunc')
    error('''fibindex_trunc'' is a required field in DTItracts')
end

endpoints   = NaN(nFib,2,3); % location of endpoints after extrapolation
endpoints_dir = NaN(nFib,2,3); % fibre direction at endpoint
ext         = NaN(nFib,2); % extension (in mm)
fibrelength = NaN(nFib,1); % fascicle length
residual    = NaN(nFib,1); % residual of polynomial fit
attach_type = NaN(nFib,2); % array indicating, for each endpoint, whether attachment is on muscle (1) or aponeurosis (2)


hwait = waitbar(0,sprintf('Extrapolating %d fibres',nFib),...
    'Name','Progress bar ExtrapolateTracts');
for fibnr = 1:1:nFib
    if any(round(linspace(1,nFib,50))==fibnr)
        waitbar(fibnr/nFib,hwait)
    end
    
    first = DTItracts.fibindex_trunc(fibnr,1);
    last  = DTItracts.fibindex_trunc(fibnr,2);
    if any(isnan([first,last]))
        PolyCoeff(fibnr,1).x  = NaN(1,order+1);
        PolyCoeff(fibnr,1).y  = NaN(1,order+1);
        PolyCoeff(fibnr,1).z  = NaN(1,order+1);        
        PolyCoeff(fibnr,1).t0 = NaN;
        PolyCoeff(fibnr,1).t1 = NaN;
        continue
    end
    if first < last
        sgn = 1;
    else
        sgn = -1;
    end
    tractpoints = [DTItracts.tracts_xyz(1,first:sgn:last)' ...
                   DTItracts.tracts_xyz(2,first:sgn:last)' ...
                   DTItracts.tracts_xyz(3,first:sgn:last)'];

    nPoints = size(tractpoints,1);
    T = zeros(nPoints,order+1);
    for o = order:-1:0
        T(:,order-o+1) = (0:nPoints-1).^o;
    end

    % Fit the polynomial
    coeff.x = (T \ tractpoints(:,1))';
    coeff.y = (T \ tractpoints(:,2))';
    coeff.z = (T \ tractpoints(:,3))';
    coeff.t0 = 0;
    coeff.t1 = nPoints-1;   
    PolyCoeff(fibnr,1) = coeff;            
    
    % Calculate residual distance of data points to polynomial fit
    if use_mex == true
        [~,dist] = FindNearestT_mex(coeff,tractpoints);
    else
        [~,dist] = FindNearestT(coeff,tractpoints);
    end
    residual(fibnr) = abs(mean(dist));
   
% Find projection of linearly extrapolated line from the endpoints to the
% surface model
        % ----- Endpoint 1 -----
        p1 =  [polyval(coeff.x,0),...
               polyval(coeff.y,0),...
               polyval(coeff.z,0)];
        d1 = [polyval(polyder(coeff.x),coeff.t0),...
              polyval(polyder(coeff.y),coeff.t0),...
              polyval(polyder(coeff.z),coeff.t0)];
        d1 = d1 / norm(d1); 

        tract_dir1 =  [polyval(coeff.x,1),...
                       polyval(coeff.y,1),...
                       polyval(coeff.z,1)] - p1;
        if sign(dot(d1,tract_dir1)) == 1
            % Change sign of direction
            d1 = -d1;
        end

        % ----- Endpoint 2 -----
        p2 =  [polyval(coeff.x,nPoints-1),...
               polyval(coeff.y,nPoints-1),...
               polyval(coeff.z,nPoints-1)];
        d2 = [polyval(polyder(coeff.x),coeff.t1),...
              polyval(polyder(coeff.y),coeff.t1),...
              polyval(polyder(coeff.z),coeff.t1)];
        d2 = d2 / norm(d2); 
        
        tract_dir2 =  [polyval(coeff.x,nPoints-2),...
                       polyval(coeff.y,nPoints-2),...
                       polyval(coeff.z,nPoints-2)] - p2;
        if sign(dot(d2,tract_dir2)) == 1
            % Change sign of direction
            d2 = -d2;
        end
        
        % Calculate intersection of projection of the endpoints with the surface model
        int_musc = intersectRaySurface(SurfModel,[p1;p2],[d1;d2],true,data_musc);
        
        if ~isempty(aponeurosis)
            % Also calculate the intersection of the projection of the
            % endpoints with the aponeurosis
            int_apo = intersectRaySurface(aponeurosis,[p1;p2],[d1;d2],false,data_apo);
            
            % Calculate distance from endpoints to intersection with muscle and
            % with the aponeurosis. Use the closest one.
            dist_to_musc = sqrt(sum((int_musc - [p1;p2]).^2,2));
            dist_to_apo  = sqrt(sum((int_apo - [p1;p2]).^2,2));
            [~,minidx] = min([dist_to_musc dist_to_apo],[],2);
            for i = 1 : 2
                if minidx(i) == 1
                    endpoints(fibnr,i,:) = int_musc(i,:);
                    attach_type(fibnr,i) = 1;
                elseif minidx(i) == 2
                    endpoints(fibnr,i,:) = int_apo(i,:);
                    attach_type(fibnr,i) = 2;
                end
            end
        else
            endpoints(fibnr,:,:) = int_musc;
            attach_type(fibnr,1:2) = 1;
        end
        endpoints_dir(fibnr,:,:) = [d1;d2];
        
%         plot3(p1(1),p1(2),p1(3),'ro','MarkerSize',8,'MarkerFaceColor','r')
%         quiver3(p1(1),p1(2),p1(3),d1(1),d1(2),d1(3),5,'r','LineWidth',2)
%         plot3(p2(1),p2(2),p2(3),'ro','MarkerSize',8,'MarkerFaceColor','c')
%         quiver3(p2(1),p2(2),p2(3),d2(1),d2(2),d2(3),5,'c','LineWidth',2)
%         
%         plot3(int_musc(:,1),int_musc(:,2),int_musc(:,3),'o','MarkerSize',8,'MarkerFaceColor','g')
%         plot3(int_apo(:,1),int_apo(:,2),int_apo(:,3),'o','MarkerSize',8,'MarkerFaceColor','m')
%         ep = min();
%         endpoints(fibnr,:,:) = ep;
        
% %%       Some code that could be uncommented for diagnostic purposes        
%         figure;
%         hold on
%         axis equal
%         patch('vertices',SurfModel.vertices,'faces',SurfModel.faces,...
%             'FaceAlpha',0.2,'FaceColor','g')
%         t_plot = linspace(coeff.t0,coeff.t1,100);
% 
%         plot3(polyval(coeff.x,t_plot),...
%             polyval(coeff.y,t_plot),...
%             polyval(coeff.z,t_plot),...
%             'LineWidth',3,'Color','b')
%         view(-34,-9)
% 
%         plot3([p1(1) endpoints(fibnr,1,1)],[p1(2) endpoints(fibnr,1,2)],[p1(3) endpoints(fibnr,1,3)],...
%             'LineWidth',2,'Color','r','MarkerSize',8,'Marker','o')
%         
%         plot3([p2(1) endpoints(fibnr,2,1)],[p2(2) endpoints(fibnr,2,2)],[p2(3) endpoints(fibnr,2,3)],...
%             'LineWidth',2,'Color','r','MarkerSize',8,'Marker','o')
% 
%%    
    % Calculate the total length of the fibre, which is the sum of the
    % extrapolated parts and the length along the polynomial-fitted curve.
    % 1) Length of polynomial fitted section (this should be very close to the
    % length of the tract on which the polynomial was fitted).
    
    [ CurveLength] = FuncCurveLength( coeff,0,nPoints-1 );
    ext(fibnr,:) = sqrt(sum(([p1;p2] - squeeze(endpoints(fibnr,:,:))).^2,2))';
    fibrelength(fibnr) = CurveLength + sum(ext(fibnr,:));

end

% Calculate the extension as percentage of the total fibre length.
pct_ext = sum(ext,2) ./ fibrelength * 100;
%% Create output arguments
% If only one output argument is provided, add outputs to the structure array
% DTItracts as fields.
if nargout == 1
    DTItracts.PolyCoeff   = PolyCoeff;
    DTItracts.fibrelength = fibrelength;
    DTItracts.ext         = ext;
    DTItracts.pct_ext     = pct_ext;
    DTItracts.endpoints   = endpoints;
    DTItracts.residual    = residual ;
    DTItracts.attach_type = attach_type ;
    DTItracts.endpoints_dir = endpoints_dir;
    varargout{1}          = DTItracts;
end
if nargout > 1
    varargout{1} = PolyCoeff;
    varargout{2} = fibrelength;
    if nargout > 2
        varargout{3} = ext;
        if nargout > 3
            varargout{4} = pct_ext;
            if nargout > 4
                varargout{5} = endpoints;
                if nargout > 5
                    varargout{6} = residual; 
                    if nargout > 6
                        varargout{7} = attach_type; 
                        if nargout > 7
                            varargout{8} = endpoints_dir; 
                        end
                    end
                end
            end
        end
    end
end
    

t_elapsed = toc;
fprintf('It took %.2f seconds to extrapolate fibre tracts.\n',t_elapsed)
close(hwait)
end % of function
