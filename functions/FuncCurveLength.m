function [ CurveLength,StraightLength ] = FuncCurveLength( P,t0,t1,varargin )
%FUNCCURVELENGTH Calculates the curve of the polynomial curve P(t) from t0 to
%t1. The curve is segmented in N pieces and the length is approximated as
%the sum of the straight-line-length of these pieces. The straight-line
%length is also calculated as the distance between the points at t0 and t1.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
% ----------------- USAGE ----------------- 
% [CurveLength, StraightLength] = CalcCurvature(P,t0,t1)
% or
% [CurveLength, StraightLength] = CalcCurvature(P,t0,t1,N)
%
% ----------------- INPUT ----------------- 
% - P: 1x1 structure array with fields x, y and z containing the
%      coefficients of the polynomial for the x, y and z-coordinates,
%      respectively.
% - t0: scalar of first t-value on the curve
% - t1: scalar of last t-value on the curve
%
% Optional:
% - N : number of segments in which the curve is subdivided to approximate
%       the curved length. If not provided, N = 1000.
%
% ----------------- OUTPUT ----------------- 
% - CurveLength: Length of the curve between t0 and t1
% - StraightLength (optional): Euclidean distance between point
%                              [x(t0) y(t0) z(t0)] and [x(t1) y(t1) z(t1)]

% Check input arguments
p = inputParser;
addRequired(p,'P',@(x) validateattributes(x,{'struct'},{'nonempty'}))
addRequired(p,'t0',@(x)validateattributes(x,{'numeric'},{'scalar'}))
addRequired(p,'t1',@(x)validateattributes(x,{'numeric'},{'scalar'}))
addOptional(p,'N',1000,@(x)validateattributes(x,{'numeric'},{'scalar'}))
parse(p, P,t0,t1,varargin{:});

N = p.Results.N;

% -- Curve length
% First, sample at a constant interval of t and calculate the segment
% lengths.
t = linspace(t0,t1,N);
dt = (t1-t0)/(N-1);
% Calculate the curve lengths
CurveLength = sum(sqrt( polyval(polyder(P.x),t).^2 +...
                        polyval(polyder(P.y),t).^2 +...
                        polyval(polyder(P.z),t).^2 ) * dt);
if nargout == 2
    StraightLength = norm([polyval(P.x,t1) polyval(P.y,t1) polyval(P.z,t1)] - ...
                          [polyval(P.x,t0) polyval(P.y,t0) polyval(P.z,t0)]);
end
end

