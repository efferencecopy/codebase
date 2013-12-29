function [varargout] = DTunpack(stro, varargin)
%
% EXAMPLE: Staircasing ==> [norms, thresh, cycperdeg, units] = DTunpack(stro)
%                 MOCS ==> [modFits, cellFits, performance] = DTunpack(stro, defaultParams)
%                QUEST ==> [thresholds, colorDirs, sfs] = DTunpack(stro, fitMeth, WINSIZE, PERFRANGE)
%
% Just a wrapper function for the unpacking functions to deal with
% staircasing and MOCS and Quest.
%
%CAH 8.20.08


if nargin < 1;
    stro = nex2stro;
end

stro = stripOutGratingTrials(stro); %remove grating catch trials!!
switch stro.sum.exptParams.expt_meth
    case 1 %MoCS
        if nargin > 1;
        	[varargout{1}, varargout{2}, varargout{3}] = DTmocsUnpack(stro, varargin{1});
        else
            [varargout{1}, varargout{2}, varargout{3}] = DTmocsUnpack(stro);
        end
    case 2 %Staircasing
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}] = DTstairsUnpack(stro);
    case 3 %QUEST
        if nargin > 2 %optional args present?
            if length(varargin) > 3
                [varargout{1}, varargout{2}, varargout{3}] = DTquestUnpack(stro, varargin{1}, varargin{2}, varargin{3});
            else
                [varargout{1}, varargout{2}, varargout{3}] = DTquestUnpack(stro, varargin{1}, varargin{2});
            end
        else
            [varargout{1}, varargout{2}, varargout{3}] = DTquestUnpack(stro, varargin{1});
        end
end

