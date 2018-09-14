function [ tn,MinDist ] = FindNearestT( P,point)
%FINDNEARESTT Finds the t-value of the point on the polynomial described in
% P that is nearest to point. Also outputs the minimal distance.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
%
O = length(P.x)-1; % Order of the polynomial
if isfield(P,'z');nDim = 3;else nDim=2;end % Dimensionality of the data

nPoints = size(point,1);
tn = NaN(1,nPoints);
MinDist = NaN(1,nPoints);
for k = 1:nPoints
    % The t-value of the minimum distance on the polynomial to 'point' can be
    % calculated by finding the minimum of the distance function of the curve
    % to 'point'. The distance function is (or actually, the squared of the
    % distance function).
    
    switch nDim
        case 2
            DistFun = conv(P.x - [zeros(1,O) point(k,1)],P.x - [zeros(1,O) point(k,1)]) + ...
                conv(P.y - [zeros(1,O) point(k,2)],P.y - [zeros(1,O) point(k,2)]);
        case 3
            DistFun = conv(P.x - [zeros(1,O) point(k,1)],P.x - [zeros(1,O) point(k,1)]) + ...
                conv(P.y - [zeros(1,O) point(k,2)],P.y - [zeros(1,O) point(k,2)]) + ...
                conv(P.z - [zeros(1,O) point(k,3)],P.z - [zeros(1,O) point(k,3)]);
    end
    
    % Now find the roots of the derivative: where the slope of the distance
    % function is zero, a minimum (or maximum) has been reached.
    t = roots(polyder(DistFun));
    
    % Only keep real part
    t = real(t);
    nT = length(t);
    
    % Sometimes multiple points are found: calculate the distance to decide
    % which one is the closest one.
    if nDim == 2
        [d,MinIndex]  = min(sum(([polyval(P.x,t) polyval(P.y,t)] - repmat(point(k,:),nT,1)).^2,2));
    elseif nDim == 3
        [d,MinIndex]  = min(sum(([polyval(P.x,t) polyval(P.y,t) polyval(P.z,t)] - repmat(point(k,:),nT,1)).^2,2));
    end
    tn(k) = t(MinIndex);
    MinDist(k) = sqrt(d);
end

end % of function

