function tf = ConnectedToViewPixx()
persistent rslt
if isempty(rslt)
    if ismac
        modelstring = 'ioreg -lw0 | grep EDID | sed "/[^<]*</s///" | xxd -p -r | strings -6';
    elseif isunix
        modelstring = 'cat /etc/X11/xorg.conf | grep -o VIEWPixx';
    else
        error('unknown operating system!');
    end
    [rc,buf] = system(modelstring);
    rslt = ~rc & strncmpi(strtrim(buf), 'viewpixx', 8);
end
tf = rslt;
