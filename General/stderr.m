function sem = stderr(data, dim)

if isempty(data)
    sem = NaN;
    return
end


if ~exist('dim', 'var')
    if ~isvector(data)
        dim = 1;
    else
        if isrow(data)
            dim = 2;
        else
            dim = 1;
        end
    end
end

sigma = std(data, [], dim);
N = size(data, dim);

sem = sigma ./ sqrt(N);
