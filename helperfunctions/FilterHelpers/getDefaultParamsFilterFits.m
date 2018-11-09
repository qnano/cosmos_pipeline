function params = getDefaultParamsFilterFits
% getDefaultParamsFilterFits get default values for ParamsFilterFits
%
% Example:
%   params = getDefaultParamsFilterFits;
%
% ParamsFilterFits: Localization Parameters
%
% Fields:
%       MaxCRLBSTD 
%       MinPValue 
%       MinPhotons
%       MinPixelDist

params.MinCRLBSTD = 0;
params.MinPValue = 0;
params.MinPhotons = 50;
params.MinBg = 0;

params.MaxCRLBSTD = 0.5;
params.MaxPValue = 1;
params.MaxPhotons = Inf;
params.MaxBg = Inf;

params.MinPixelDist = 0;

end

