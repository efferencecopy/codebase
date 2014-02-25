function [stack info_tmp] = colorMerge(params)

% EXAMPLE INPUT PARAMS
%
%   params.mouse = 'CH_112613_B';
%   parms.objective = '2x';
%   params.contrastMethod = 'none';
%   params.npix = 0;

% cd to where the images are
global GL_DATPATH
cd([GL_DATPATH, filesep, params.mouse, filesep, 'Histology', filesep, 'Raw Images']);

% grab the names in the directory
d = dir;

% load in (and update) the mouseDB
mdb = initMouseDB();

% initialize the structure of images
[img.green, img.red, img.blue] = deal({});


for a = 1:numel(d);
    
    % display the progress
    if a == 1 || ~rem(a,5)
        fprintf('%d more images to unpack\n', numel(d)-(a-1));
    end
    
    if ~any(regexpi(d(a).name, '^[\.]|thumbs')) % skip the hidden files

        
        % make sure the objective used is correct
        sliceObjective = regexpi(d(a).name , '_\d+x', 'match');
        sliceObjective = sliceObjective{1}(2:end);
        if strcmpi(sliceObjective, params.objective)
            
            
            % id the plate number
            plate = regexpi(d(a).name , '_p\d+', 'match');
            plate = str2double(plate{1}(3:end));
            
            % id the slice number
            slice = regexpi(d(a).name , '_s\d+', 'match');
            slice = str2double(slice{1}(3:end));
            
            
            % red? or green?
            if regexp(d(a).name, '_red'); color = 'red'; end
            if regexp(d(a).name, '_green'); color = 'green'; end
            if regexp(d(a).name, '_blue');   color = 'blue'; end
            if regexp(d(a).name, '_white');   continue; end
            
            
            % unpack the images
            img_tmp = imread(d(a).name);
            info_tmp = imfinfo(d(a).name);
            switch info_tmp.ColorType
                case 'truecolor'
                    img_tmp = preProcessImg(img_tmp, info_tmp, params.npix, params.contrastMethod);
                case 'grayscale'
                    % no need to do anything, already in grayscale
            end
            
            
            % put the image in a structure according to it's position in
            % the brain, and the color channel
            img = setfield(img, color, {plate, slice},  {img_tmp});
            
        end
    end
end



% figure out some histology parameters for this mouse (used below):
mdbidx = regexpi({mdb.mice(:).name}', params.mouse);
mdbidx = ~cellfun(@isempty, mdbidx);
thickness = mdb.mice(mdbidx).histo.thickness;
slicesPerPlate = mdb.mice(mdbidx).histo.slicesPerPlate;


% combine the images into a merged truecolor RGB. Arrange them in a stack.
idx = 1;
numPlates = max([size(img.red,1), size(img.green,1), size(img.blue,1)]);
numSlices = max([size(img.red,2), size(img.green,2), size(img.blue,2)]);

for p = 1:numPlates
    for sl = 1:numSlices
        
        % find the appropriate images
        if ((size(img.red,1)>=p) && (size(img.red,2)>=sl)) && ~isempty(img.red{p,sl})
            ch_red = img.red{p, sl};
        else
            ch_red = zeros(size(img_tmp));
        end
        
        if (size(img.green,1)>=p) && (size(img.green,2)>=sl) && ~isempty(img.green{p,sl})
            ch_green = img.green{p, sl};
        else
            ch_green = zeros(size(img_tmp));
        end
        
        if (size(img.blue,1)>=p) && (size(img.blue,2)>=sl) && ~isempty(img.blue{p,sl})
            ch_blue = img.blue{p, sl};
        else
            ch_blue = zeros(size(img_tmp));
        end
        
        % do some basic error checking
        if isempty(ch_green) && ~isempty(ch_green);
            error('One channel is defined but not the other')
        end
        if all([ch_green(:);ch_red(:);ch_blue(:)]==0); continue; end % all channels are undefined... no big deal so no error
        
        
        % add images to the stack, make up a fake blue channel. During
        % image acquisition, the microscope objective flips the image L/R
        % and U/D, so flip all of them back...
        merge = cat(3, ch_red, ch_green, ch_blue);
        for a = 1:3
            merge(:,:,a) = rot90(merge(:,:,a), 2);
        end
        stack.img{idx} = merge;
        
        % figure out where the slice was in the brain based off of the
        % slice thickness and the slice number
        stack.loc(idx) = sum(slicesPerPlate(1:p-1).*thickness) + ((sl-1).*thickness);

        % update the index
        idx = idx + 1;
        
    end
end



    

    