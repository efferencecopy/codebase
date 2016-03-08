function filename = generate_filename(k, tag)
% we try to contain as much data as possible in the filename
% in order for it to be relly unique

    X_min = min(k.X(:)); Y_min = min(k.Y(:));
    X_max = max(k.X(:)); Y_max = min(k.Y(:));
    gdX = k.X(2, 1) - k.X(1, 1);
    gdY = k.Y(2, 1) - k.Y(1, 1);

    X_src_min = min(k.X_src(:)); Y_src_min = min(k.Y_src(:));
    X_src_max = max(k.X_src(:)); Y_src_max = max(k.Y_src(:));
    gdX_src = k.X_src(2, 1) - k.X_src(1, 1);
    gdY_src = k.Y_src(1, 2) - k.Y_src(1, 2);
    
    X_el_pos_min = min(k.el_pos(:,2));
    X_el_pos_max = max(k.el_pos(:,2));
    Y_el_pos_min = min(k.el_pos(:,1));
    Y_el_pos_max = max(k.el_pos(:,1));
    X_el_pos_mean = mean(k.el_pos(:,2));
    Y_el_pos_mean = mean(k.el_pos(:,1));
    
    %num2str(), ' '

    filename = [num2str(k.R), '_', num2str(k.h), '_',...
        num2str(k.sigma), '_', num2str(X_min), '_', num2str(X_max), '_',...
        num2str(Y_min), '_', num2str(Y_max), '_', num2str(gdX), '_',...
        num2str(gdY), '_', num2str(X_src_min), '_', num2str(X_src_max), '_',...
        num2str(Y_src_min), '_', num2str(Y_src_max), '_', num2str(gdX_src),...
        '_', num2str(gdY_src), '_', num2str(X_el_pos_min),'_', num2str(X_el_pos_max),'_',...
        num2str(Y_el_pos_min), '_', num2str(Y_el_pos_max),...
        num2str(X_el_pos_mean), num2str(Y_el_pos_mean),...
        length(k.X(:)), length(k.Y(:)), length(k.el_pos(:,1))];
    filename(filename=='.') = 'x';
    s = which('kcsd2d.m');
    s(length(s)-8:length(s)) = [];
    filename = [s,'/data/', tag, '_', filename];
end
    
