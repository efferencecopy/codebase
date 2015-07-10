function sem = sem(data, dim)

if ~isvector(data)
    if ~exist('dim', 'var')
        dim = 1;
    end
else
    if isrow(data)
        dim = 2;
    else
        dim = 1;
    end
end

sigma = std(data, [], dim);
N = size(data, dim);

sem = sigma ./ sqrt(N);
