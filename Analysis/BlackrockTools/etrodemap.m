function map_of_blk_ch_locations = etrodemap(etrodetype)


switch etrodetype
    case 'a32-4x8'
        
        load etrodeMapping_A32_4x8.mat
        firstShank = [5;4;6;3;7;2;8;1];
        mtx = bsxfun(@times, ones(8,4), [0,8,16,24]);
        map_of_contact_locations = bsxfun(@plus, mtx, firstShank);
        map_of_blk_ch_locations = nan(size(map_of_contact_locations));
        for i_site = 1:numel(map_of_contact_locations)
            ch = map_of_contact_locations(i_site);
            idx = A32_4x8.electrodeNumber == ch;
            blk_channel = A32_4x8.headstageChannel(idx);
            map_of_blk_ch_locations(i_site) = blk_channel;
        end
            
    otherwise
        error('Electrode type not recognized')
end
