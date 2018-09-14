function [h,h_norm] = plot_surface( surf_model,varargin )
%PLOT_SURFACE Plots a surface model with the default settings, or custom
%settings when provided.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% March 2018
% 
% Usage:
% handle = plot_surface(surf_model,'parameter',<value>)
% 
% surf_model : structure containing fields vertices and faces or filename 
%              of an STL-file.
% Optional inputs:
% FaceColor  : color of the faces of the surface model. Default: 'r'
% EdgeColor  : color of the edge of the surface model. Default: 'none'
% FaceAlpha  : transparency value of the faces. Default: 0.3
% EdgeAlpha  : transparency value of the edges. Default: 1
% ShowNormals : Show normal vectors of the surface model.


p = inputParser;
addRequired(p,'surf_model')
addParameter(p,'FaceColor','r')
addParameter(p,'FaceAlpha',0.3,@isscalar)
addParameter(p,'EdgeColor','none')
addParameter(p,'EdgeAlpha',1,@isscalar)
addParameter(p,'ShowNormals',false)

parse(p,surf_model,varargin{:})

% Read the STL-file 
if ~isstruct(surf_model)
    surf_model = stlread(surf_model);
end

holdstate = ishold(gca);
hold on
h = patch('Vertices',surf_model.vertices,...
               'Faces',surf_model.faces,...
               'FaceColor',p.Results.FaceColor,...
               'FaceAlpha',p.Results.FaceAlpha,...
               'EdgeColor',p.Results.EdgeColor,...
               'EdgeAlpha',p.Results.EdgeAlpha);

if p.Results.ShowNormals == true
    if ~isfield(surf_model,'normals')
        surf_model.normals = facenormals(surf_model);
    end
    
    C = squeeze(mean(reshape(surf_model.vertices(surf_model.faces,:)',3,size(surf_model.faces,1),3),3))';
    h_norm = quiver3(C(:,1),C(:,2),C(:,3),...
        surf_model.normals(:,1),...
        surf_model.normals(:,2),...
        surf_model.normals(:,3),0,...
        'Color','r',...
        'AutoScale','off');
end

axis equal
if holdstate == 0
    hold(gca,'off')
else
    hold(gca,'on')
end

end

