function sem = stderr(data, dim)

if isempty(data)
    sem = NaN;
    return
end


if ~exist('dim', 'var')
    if isvector(data)
        if isrow(data)
            dim = 2;
        else
            dim = 1;
        end
    else
        dim = 1;
    end
end

sigma = std(data, [], dim);
N = size(data, dim);

sem = sigma ./ sqrt(N);
