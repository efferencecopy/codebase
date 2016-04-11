function f = myfig(position, name)

f = figure;
if exist('position', 'var')
    set(f, 'position', position')
end
if exist('name', 'var')
    set(f, 'name', name)
end