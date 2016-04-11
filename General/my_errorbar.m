function [h_line, h_errbar] = my_errorbar(xx, yy, errUp, varargin)

if ~exist('errDown', 'var')
    errDown = errUp;
end

assert(isvector(errUp), 'ERROR, my_errbar only works with vector inputs')

% deal with the optional lineseries inputs
if ~isempty(varargin)
    linespec = [];
    for i_spec = 1:numel(varargin)
        
        tmp = varargin{i_spec};
        if isnumeric(tmp)
            if isscalar(tmp)
                tmp = num2str(tmp);
            elseif isvector(tmp)
                tmp = sprintf('[%s]',num2str(tmp));
            end
        elseif ischar(tmp)
            tmp = sprintf('''%s''', tmp);
        end
        linespec = [linespec, ',' tmp];
    end
end


% plot the line series
hold on,
plotstring = 'h_line = plot(xx, yy';
if isempty(varargin)
    plotstring = [plotstring, ')'];
else
    plotstring = [plotstring, linespec, ');'];
end
eval(plotstring);


% plot the error bars
x_err = [xx(:)'; xx(:)'];
y_err = [yy(:)' + errUp(:)' ; yy(:)' - errDown(:)'];
plotstring = 'h_errbar = plot(x_err, y_err';
if isempty(varargin)
    plotstring = [plotstring, ',''marker'', ''none'');'];
else
    plotstring = [plotstring, linespec, ',''marker'', ''none'');'];
end
eval(plotstring);


