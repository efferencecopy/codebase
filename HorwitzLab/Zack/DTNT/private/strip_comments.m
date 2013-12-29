function str_out = strip_comments(str_in)
str_out = '';
if ~isempty(str_in)
    percent_idx = strfind(str_in, '%');
    if ~isempty(percent_idx)
        str_out = str_in(1:percent_idx(1)-1);
    else
        str_out = str_in;
    end
end
