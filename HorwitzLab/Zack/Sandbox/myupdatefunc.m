function txt = myupdatefunc(~,event_obj)
zout = evalin('base', 'zout;');
pos = get(event_obj,'Position');
sse = zout{3}(pos(2),pos(1));
params = zout{2}{pos(2),pos(1)};

txt = {sprintf('[i j] = [%d %d]', pos(2), pos(1)), ...
        sprintf('sse: %f', sse), ...
        sprintf('params:\n%f\n%f\n%f', params)};