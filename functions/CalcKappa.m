function kappa = CalcKappa( P,t )
%CALCKAPPA Calculates the curvature kappa of the 3D polynomial curve P at
%points [x(t),y(t),z(t)] using the Frenet formula.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
% ----------------- USAGE ----------------- 
% kappa = FuncCurvature(P,t)
% 
% ----------------- INPUT ----------------- 
% - P: nx1 structure array with fields x, y and z containing the
%      coefficients of the polynomial for the x, y and z-coordinates,
%      respectively.
%
% ----------------- OUTPUT ----------------- 
% - kappa: n x 1 array of curvature values

kappa=zeros(1,length(t));
for i = 1:length(t)
    tc = t(i);
    % First derivative to t
    rdot = [ polyval(polyder(P.x),tc);...
             polyval(polyder(P.y),tc);...
             polyval(polyder(P.z),tc)];

    % Second derivative to t
    rddot = [polyval(polyder(polyder(P.x)),tc);...
             polyval(polyder(polyder(P.y)),tc);...
             polyval(polyder(polyder(P.z)),tc)];
         
% The curvature is defined as the norm of the cross product of the first
% and second derivative of P, divided by the norm of the first derivative
% to the power 3 (Frenet formula).
    kappa(i) = norm(cross(rdot,rddot)) / norm(rdot)^3;
end

end % of function

