%% Data files for EB_012214_B cell #1
fin


fnamePrefix = '2014_02_14_';
fnameSuffix = 2:1:24;
 % group the files according to whether or not the AS, FS, iris was open.
 % 1= All open. 2=FS closed, 3=FS+AS closed
 group = [2,2,2,1,2,2,1,1,1,2,2,2,2,1,2,1,2,3,3,2,2,3,2];
 distance = [0, 200, 300, 300, 400, 500, 500, 750, 1000, 250, 300, 350, 400, 400, 450, 450, 300, 300, 250, 250, 200, 200, 200];
 
 % analysis windows
 plt_window = [1.035 1.14];
 Im_window = [1.045 1.048];
 base_window = [1.03 1.037];
 led_window = [1.04 1.043];
 
 
%% Data files for EB_012214_B cell #2
fin


fnamePrefix = '2014_02_14_';
fnameSuffix = 28:1:32;

 % group the files according to whether or not the AS, FS, iris was open.
 % 1= All open. 2=FS closed, 3=FS+AS closed
 group = [1,1,2,1,3];
 distance = [0,200,200,200,200];
 
  % analysis windows
 plt_window = [1.035 1.14];
 Im_window = [1.041 1.042];
 base_window = [1.03 1.037];
 led_window = [1.04 1.043];
 
 %% Data files for EB_012414_A cell #5
fin


fnamePrefix = '2014_02_13_';
fnameSuffix = [17:1:23,25];

 % group the files according to whether or not the AS, FS, iris was open.
 % 1= All open. 2=FS closed, 3=FS+AS closed
 group = [1,3,3,3,1,1,3,3];
 distance = [0,0,0,150,150,300,300,300];
 
 % analysis windows
 plt_window = [1.33 1.42];
 Im_window = [1.34 1.36];
 base_window = [1.26 1.33];
 led_window = [1.34 1.343];
 
 
 

 
 %% ANALYZE AND PLOT

 for a = 1:numel(fnameSuffix);
     
     num = fnameSuffix(a);
     suffix = [repmat('0', 1, 3-double(num>9)), num2str(num)];
     fname = [fnamePrefix, suffix];
     
     fpath = findfile(fname, GL_DATPATH, '.abf');
     [dat, si, h] = abfload(fpath);
     
     N = size(dat, 1);
     tt = (0:N-1) .* (si * 1e-6);
     
     % plot the data
     window = plt_window;
     idx = (tt >= window(1)) & (tt < window(2));
     figure
     set(gcf, 'position', [440     9   783   808])
     subplot(2,1,2)
     set(gca, 'OuterPosition', [-2.77556e-17 0.0875229 1 0.204338])
     plot(tt(idx), permute(dat(idx,2,:), [1,3,2]))
     axis tight
     ylabel('voltage')
     xlabel('time (sec)')
     subplot(2,1,1)
     plot(tt(idx), permute(dat(idx,1,:), [1,3,2]))
     ylabel('Current')
     axis tight
     set(gca, 'OuterPosition', [-2.77556e-17 0.300049 1 0.675623])
     title(sprintf(fname))
     
     
     % deal with the Im
     window = Im_window; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     Im_min{a} = mean(tmp, 1);
     
     % deal with the baseline
     window = base_window; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     baseline = mean(tmp, 1);
     Im_min{a} = Im_min{a} - baseline;
     
     % deal with the LED voltage
     window = led_window; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,2,:), [1,3,2]);
     V_led{a} = mean(tmp, 1);  
     
 end
     
 
 
 % plot each of the groups
 for i = 1:3
     l_group = group == i;
     if sum(l_group)==0; continue; end
     dist_grp = distance(l_group);
     
     clr = colormap(jet);
     idx = round(linspace(1, size(clr,1), sum(l_group)));
     clr = clr(idx,:);
     
     figure, hold on,
     idx = find(l_group);
     for a = 1:sum(l_group);
         plot(V_led{idx(a)}, Im_min{idx(a)}, 'color', clr(a,:), 'linewidth', 3)
     end
     xlabel('voltage')
     ylabel('evoked current')
     
     switch i
         case 1
             title('Open Iris')
         case 2
             title('FS closed')
         case 3
             title('FS, AS closed')
     end
    
     % deal with the legend
     distAsCell =  mat2cell(distance(l_group)', ones(1,sum(l_group)))';
     distForLeg = cellfun(@num2str, distAsCell, 'uniformoutput', false);
     legend(distForLeg)
     
 end
 
 