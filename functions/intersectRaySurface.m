function I = intersectRaySurface(FV,pts,D,inside,opt)
%%INTERSECTRAYSURFACE calculates the location of a point projected on the
%%surface FV (with fields 'faces' and 'vertices') along direction d.
%
% Bart Bolsterlee, Neuroscience Research Australia
% February 2017
%
% ----------------- INPUT ----------------- 
% FV     : structure array with fields 'faces' and 'vertices'
% pts    : n x 3 matrix with xyz coordinates of points
% D      : n x 3 matrix of direction along which the points will be projected
% 
% Optional 4th input argument:
% inside: logical (true/false, 0/1) indicating whether the point is
% initially inside (true/1) or outside (false/0) the surface model. If
% inside = true/false, the projection of points outside/inside the surface 
% model is not calculated. The initial point is then returned. Default:
% true
%
% Optional 5th input argument:
% opt  : structure containing fields V0,V1,V2, U,V, N, UU,UV,VV and DD with
%        edge vectors/dot products of edge vectors and normals (see below).
%
% ----------------- OUTPUT ----------------- 
% I = n x 3 matrix with xyz coordinates of projections of the points on the
% surface
%
% The algorithm to calculate the intersection of a ray and a triangle was
% inspired by some linear algebra wizardry 'stolen' from the c-code by
% Dan Sunday available on:
% http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()

if nargin <= 3
    if isempty(inside)
        inside = true;
    end
end

if nargin <= 4
    % Get all vertices locations
    opt.V0 = FV.vertices(FV.faces(:,1),:);
    opt.V1 = FV.vertices(FV.faces(:,2),:);
    opt.V2 = FV.vertices(FV.faces(:,3),:);

    % Calculate edge vectors
    opt.U = opt.V1 - opt.V0; 
    opt.V = opt.V2 - opt.V0;

    % Calculate dot product of edge vector
    opt.UU = sum(opt.U.*opt.U,2);
    opt.UV = sum(opt.U.*opt.V,2);
    opt.VV = sum(opt.V.*opt.V,2);
    opt.DD = opt.UV.*opt.UV - opt.UU.*opt.VV;

    % Calculate normal of faces
    opt.N = cross(opt.U,opt.V);
end

IN_SURF = inside_surface(FV,pts,opt);

% Calculate which points are inside the surface
nF = size(FV.faces,1);
I  = zeros(size(pts,1),3);
for i = 1 : size(pts,1)
    point = pts(i,:);
%     in_surf = inside_surface(FV,point);
    
    if (~IN_SURF(i) && inside) || (IN_SURF(i) && ~inside)
        % If the point is already outside/inside the surface (depending on 
        % the value of 'inside'), keep the original point as the projection.
        I(i,1:3) = point;
        continue
    end
    % Vector implementation of the code available at:
    % http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
    W0 = ones(nF,1) * point - opt.V0;
    A = -sum(opt.N .* W0,2);
    B =  sum(opt.N .* (ones(nF,1) * D(i,:)),2);
    
    % if R < 0 the ray points away from the triangle and no intersection
    % will be found. Exclude all triangles the ray points away from.
    R = A ./ B;
    included = R >= 0;
    nIncl = sum(included);
    
    % Calculate intersection of ray and triangle planes
    II = ones(nIncl,1) * point + (R(included)*ones(1,3)) .* (ones(nIncl,1) * D(i,:));
    
    % Calculate which triangles are intersected
    W = II - opt.V0(included,:);
    WV = sum(W.*opt.V(included,:),2);
    WU = sum(W.*opt.U(included,:),2);
    S = (opt.UV(included) .* WV - opt.VV(included) .* WU) ./ opt.DD(included);
    T = (opt.UV(included) .* WU - opt.UU(included) .* WV) ./ opt.DD(included);
    
    idx = find(S >= 0 & S <= 1 & T >= 0 & (S+T) <= 1);
    if length(idx) == 1
        I(i,1:3) = II(idx(1),1:3); % only one intersection found
    elseif isempty(idx)
        I(i,1:3) = NaN(1,3); % no intersection found
    elseif length(idx) > 1
        % multiple intersections found. Select the closest one.
        distance = sqrt(sum((II(idx,1:3) - ones(length(idx),1) * point).^2,2));
        [~,closest_idx] = min(distance);
        I(i,1:3) = II(idx(closest_idx),1:3);
    end
        
%     %% Some code for diagnostic purposes
%     figure;
%     hold on
%     axis equal
%     patch('vertices',FV.vertices,'faces',FV.faces,...
%         'FaceAlpha',0.2,'FaceColor','g')
%     plot3(II(idx,1),II(idx,2),II(idx,3),'ro',...
%         'MarkerFaceColor','r','MarkerEdgeColor','none',... %
%     'MarkerSize',8)
%     quiver3(point(1),point(2),point(3),...
%         D(i,1),D(i,2),D(i,3),10,'-','filled','LineWidth',3,...
%         'Color','r')
%     plot3([point(1) II(idx,1)],[point(2) II(idx,2)],[point(3) II(idx,3)],...
%         'LineWidth',2,'Color','b')
% 
%     view(-34,-9)
% %     

end % of function




