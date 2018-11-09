function params = getDefaultParamsFit
% getDefaultParamsFit get default values for ParamsFit
%
% Example:
%   params = getDefaultParamsFit;
%
% ParamsFit: Localization Parameters 
%
% Fields:
%       FitSigma
%       BoxSize
%       Iterations 
%       MaxCudaFits
%       PSFSigma

 params.FitSigma = false;
 params.Iterations = 10;
 params.MaxCudaFits = 100000;
 params.PSFSigma=1.39;
 params.BoxSize=5;%3*(2*params.PSFSigma+1);
 end

