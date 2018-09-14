function curvature = CalcCurvature( PolyCoeff)
%CALCCURVATURE calculates the curvature of all fibres in input 
% structure PolyCoeff. The curvature is the mean curvature of 100 points
% sampled along the polynomial curve at equal intervals of t.
%
% Bart Bolsterlee, Neuroscience Research Australia (NeuRA)
% February 2017
% ----------------- USAGE ----------------- 
% kappa = CalcCurvature(PolyCoeff)
%
% ----------------- INPUT -----------------
% - PolyCoeff : n x 1 structure array with fields x,y, z, t0 and t1 (output
%               of ExtrapolateTracts)
%         
% ----------------- OUTPUT -----------------
% - curvature : n x 1 array of mean curvatures per fibre.
%
nFib = length(PolyCoeff);
hwait = waitbar(0,sprintf('Calculating curvature for %d fibres',nFib),...
    'Name','Progress bar CalcCurvature');

persistent ShowWarningMex
if isempty(ShowWarningMex);ShowWarningMex = true;end
if exist('CalcKappa_mex','file') == 3
    use_mex = true;
else
    use_mex = false;
    if ShowWarningMex == true
        warning('The MEX-file CalcKappa_mex is not available on the MATLAB path. Curvature calculation may be slow.')
        ShowWarningMex = false; % only show the warning message once
    end
end
tic
curvature = NaN(nFib,1);
for fibnr = 1:nFib
    if any(round(linspace(1,nFib,11))==fibnr)
        waitbar(fibnr/nFib,hwait)
    end
    
    if isempty(PolyCoeff(fibnr).x)
        continue
    end

    % sample at 100 equal intervals of t
    t = linspace(PolyCoeff(fibnr).t0,...
        PolyCoeff(fibnr).t1,100);
    if use_mex == true
        curv = CalcKappa_mex(PolyCoeff(fibnr),t');
    else
        curv = CalcKappa(PolyCoeff(fibnr),t');
    end
    curvature(fibnr) = mean(curv) * 1000; % convert from 1/mm to 1/m
end
t_elapsed = toc;
fprintf('It took %.2f seconds to calculate curvatures.\n',t_elapsed)

close(hwait)
end

