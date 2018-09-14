function IN = inside_surface(FV,pts,opt)
%%INSIDE_SURFACE calculates whether the points given in pts (n x 3 matrix
%%of x y z coordinates) are inside (1) or outside (0) the triangulated
%%surface given by FV (containing fields faces and vertices).
%
% Optionally, the edge vectors and dot products of the edge vectors can be
% provided as an additional input structure 'opt'. This speeds up
% calculation when inside_surface is used in a loop.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
% This algorithm projects a 'ray' from the query point in a given direction
% and counts how many triangles of the surface model are intersected by the
% ray. An odd number of intersections means the point is inside the surface
% and an even number means the point is outside the surface.
%
% The algorithm to calculate the intersection of a ray and a triangle was
% inspired by some linear algebra wizardry which I stole from the c-code by
% Dan Sunday available on:
% http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()

if nargin < 3
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

nP = size(pts,1);
NINT = zeros(nP,1); % Number of intersections
for i = 1 : size(pts,1)
    point = pts(i,:);
    
    % The ray is projected along the x-axis from point [px py pz].
    % Triangles that are completely below or above the plane y = py and z =
    % pz can never intersect with this ray, so these faces can be excluded
    % now. This eliminates most faces from the surface and will make the
    % computation much faster.
    
    d  = [1 0 0]; % ray direction vector in positive x-direction
    excl = (opt.V0(:,2) < pts(i,2) & opt.V1(:,2) < pts(i,2) & opt.V2(:,2) < pts(i,2)) | ... % all vertices below y = yp (i.e. triangle is completely below plane y = yp)
           (opt.V0(:,2) > pts(i,2) & opt.V1(:,2) > pts(i,2) & opt.V2(:,2) > pts(i,2)) | ... % all vertices above y = yp
           (opt.V0(:,3) < pts(i,3) & opt.V1(:,3) < pts(i,3) & opt.V2(:,3) < pts(i,3)) | ... % all vertices below z = zp
           (opt.V0(:,3) > pts(i,3) & opt.V1(:,3) > pts(i,3) & opt.V2(:,3) > pts(i,3));      % all vertices above z = zp
    %
    %     figure;hold on axis equal patch('vertices',FV.vertices,'faces',FV.faces,...
    %         'FaceAlpha',0.2,'FaceColor','g')
    %     patch('vertices',FV.vertices,'faces',FV.faces(~excl,:),...
    %         'FaceAlpha',1,'FaceColor','r')
    %     %     plot3( pts(i,1), pts(i,2), pts(i,3),'ro',... %
    %     'MarkerFaceColor','r','MarkerEdgeColor','none',... %
    %     'MarkerSize',8) quiver3( pts(i,1), pts(i,2), pts(i,3),...
    %         d(1),d(2),d(3),10,'-','filled','LineWidth',3,... 'Color','r')
    %
    %     view(-34,-9)
    
    % For the leftover faces, calculate whether the ray intersects the
    % face. The following code is the vector implementation of the c-code
    % given at.:
    % http://geomalgorithms.com/a06-_intersect-2.html#intersect3D_RayTriangle()
    nF = sum(~excl);
    
    W0 = ones(nF,1) * point - opt.V0(~excl,:);
    A = -sum(opt.N(~excl,:) .* W0,2);
    B =  sum(opt.N(~excl,:) .* (ones(nF,1) * d),2);
    
    % if R < 0 the ray points away from the triangle and no intersection
    % will be found.
    R = A ./ B;
    
    % calculate intersection of ray and plane
    II = ones(nF,1) * point + (R*ones(1,3)) .* (ones(nF,1) * d);
    
    W = II - opt.V0(~excl,:);
    WV = sum(W.*opt.V(~excl,:),2);
    WU = sum(W.*opt.U(~excl,:),2);
    S = (opt.UV(~excl,:) .* WV - opt.VV(~excl,:) .* WU) ./ opt.DD(~excl,:);
    T = (opt.UV(~excl,:) .* WU - opt.UU(~excl,:) .* WV) ./ opt.DD(~excl,:);
    
    % Count how many times the surface was intersected
    NINT(i) = sum((S >= 0 & S <= 1 & T >= 0 & (S+T) <= 1) & R >= 0);
    
end

% An odd/even number of intersections means the points is inside/outside.
IN = mod(NINT,2)~=0;

end % of function




