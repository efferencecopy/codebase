function lines = open_parse_filenames(textfile)
fid = fopen(textfile);
raw_lines = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);

raw_lines = raw_lines{1};
lines = cellfun(@strip_comments, raw_lines, 'unif', 0);
lines = strtrim(lines(~cellfun(@isempty, lines)));
