function [ normals ] = facenormals( FV )
%FACENORMALS  normal vectors to the surface defined by the
%structure FV containing vertices and faces.

% Get all edge vectors
e1=FV.vertices(FV.faces(:,1),:)-FV.vertices(FV.faces(:,2),:);
e3=FV.vertices(FV.faces(:,3),:)-FV.vertices(FV.faces(:,1),:);
% 
% Calculate normal of face
 normals = cross(e3,e1);
 % Normalize to unit vector
 normals = normals ./ repmat(sqrt(normals(:,1).^2+normals(:,2).^2+normals(:,3).^2),1,3); 
end

