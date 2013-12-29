function [chanOrder, relDist] = reorderPinDistance(POSITION, BANK)

anchorOrderMapping = {'ul', 'ur', 'll', 'lr', 'up', 'right', 'center'};
possibleAnchors    =  [91,   100,  1,    10,   95,   60,     55 ];
anchor = possibleAnchors(strcmpi(anchorOrderMapping, POSITION));

array_layout = flipud(reshape(1:100,[10 10])');
pinPositions = pins_rel2abs(BANK, 1:32);

relDist = nan(numel(pinPositions), 1);
for a = 1:numel(pinPositions)
    [r, c] = find((array_layout == anchor) | (array_layout == pinPositions(a)));
    relDist(a) = sqrt( (abs(r(1)-r(2))^2) + (abs(c(1)-c(2))^2) );
end

[~, chanOrder] = sort(relDist);