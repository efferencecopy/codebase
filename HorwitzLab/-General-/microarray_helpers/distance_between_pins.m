% Pass in absolute-numbered pins and get the pair-wise distance between all
% possible pairs. The output is a column vector whose n^2 entries are:
%  [d(p1,p1) d(p1,p2) ... d(p1,pn) d(p2,p1) ... d(p2,pn) ... d(pn,pn)]'
%  where p1,...,pn are your inputs to this function (i.e., preserves your
%  order while taking all pairs). You can also pass in two pins as two
%  arguments to this function and get the distance between them.
function dist = distance_between_pins(varargin)
persistent array_layout
array_layout = flipud(reshape(1:100,[10 10])');
dist = [];
if nargin == 2
    p1 = varargin{1}; p2 = varargin{2};
    [ii,jj] = find(array_layout==p1 | array_layout==p2);
    dist = sqrt((ii(2)-ii(1))^2 + (jj(2)-jj(1))^2);
elseif nargin == 1
    [~,loc] = ismember(array_layout, varargin{1});
    [~,orig_order] = sort(loc(loc ~= 0));
    [ii,jj] = find(loc); gridpoints = [ii jj];
    gridpoints = gridpoints(orig_order,:); % preserve the user's order
    [n,d] = size(gridpoints);
    dist = sqrt(sum(bsxfun(@minus,reshape(gridpoints,[n 1 d]),reshape(gridpoints,[1 n d])).^2,3));
    dist = dist(:);
end
