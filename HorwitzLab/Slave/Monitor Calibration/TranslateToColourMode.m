function out = TranslateToColourMode(in, oldstyle16bit, ccmode)
% out = TranslateToColourMode(in [, oldstyle16bit] [, ccmode])
%
% Takes an NxMxP (P=1,3) matrix and converts it to an Nx2MxP matrix that is
% appropriate for driving the Bits++/ViewPixx in colour mode. In this mode, the
% Bits++/ViewPixx expects the image matrix to have a certain structure, and this
% function ensures that structure.
%
% This function is wise to the PsychImaging pipeline, which automates some of
% the nuts and bolts of driving a high-resolution (color bit depth) display
% device. In particular, it will poll the color conversion mode from
% Psychtoolbox's Bits++ driver. This `ccmode` determines how to modify the
% inputs to drawing commands (Fill*, (Make|Draw)Texture, etc.) to achieve the
% highest color bit depth possible. The color conversion mode is only used when
% the user invokes the PsychImaging pipeline.
%
% For backwards compatibility, this function accepts image matrices with
% intensities encoded as DAC values (in the range 0-255 or 0-65535). Assuming
% PsychImaging wasn't used, `oldstyle16bit` of 1 encodes the image matrix by
% splitting the 16-bit number into adjacent high and low bytes (which the Bits++
% uses to achieve 14-bit color), or (flag set to 0) puts the 8-bit number into
% the high byte and zeros the low bytes. This flag is forced to 1 if a single
% image element has intensity > 255.
%
% This function can accept image matricies as uint8s, uint16s, or doubles and
% returns a matrix of the same class.
%
% If the PsychImaging pipeline or the `ccmode` argument was used, this code will
% process the image (intensities 0.0-1.0) based on the current window's
% specified color conversion mode. This entails modifying the image matrix by
% adding blank columns (mode 1), repeating every column (mode 2), or doing
% nothing (mode 0). You can read more about these modes through `help
% PsychImaging` and looking for "EnableBits++Color++Output."
%
% With PsychImaging, set the 6th argument `floatprecision` to 2 when making
% textures:
%    texture = Screen('MakeTexture', win, TranslateToColourMode(im), [], [], 2).
%
% GDLH 2007/09/05 - Made the initial version.
% ZALB 2013/04/16 - Support for Bits++ conversion modes in PsychImaging.

persistent warn_user
warn_user = true;

if ndims(in) < 2 || ndims(in) > 3 || ~ismember(size(in,3), [1 3])
    error('Input argument has to be NxMxP (P=1,3)');
end

if nargin < 2 || isempty(oldstyle16bit) || ~ismember(oldstyle16bit, [0 1])
    oldstyle16bit = 0;
end

if nargin < 3 || isempty(ccmode)
    ccmode = -1;
end

if ccmode == -1
    open_windows = Screen('Windows');
    if ~isempty(open_windows)
        try
            ccmode = BitsPlusPlus('GetColorConversionMode', max(open_windows));
        catch %#ok<CTCH>
            if warn_user
                warn_user = false;
                warning(['Failed to query the current color conversion mode. I''m assuming you ' ...
                    'aren''t using the PsychImaging pipeline to draw high-resolution textures.']); %#ok<WNTAG>
            end
        end
    end
end

out = TranslateToColourModeMex(in, oldstyle16bit, ccmode);
