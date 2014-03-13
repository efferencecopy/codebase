%% DATA FILES FOR REVERSAL POTENTIAL ANALYSIS.

fin

mdb = initMouseDB();

in = {'EB_012414_A', '2014_02_13_0002', -60, 1.3306, 1.3310;...
      'EB_012414_A', '2014_02_13_0014', -60, 1.3296, 1.3302;...
      'EB_012414_A', '2014_02_13_0029', -60, 1.3320, 1.3331;...
      'EB_012414_A', '2014_02_13_0038', -60, 1.3294, 1.3304;...
      'EB_012414_A', '2014_02_13_0039', -70, 1.3302, 1.3319;...
      'EB_012214_C', '2014_02_11_0009', -60, 1.3333, 1.3359;...
      'EB_012214_A', '2014_02_09_0002', -60, 1.3313, 1.3333};
  
  
% roll through the data files, plot the raw data, and extract the
% parameters that are important

for a = 1:size(in,1)
    
    
    fpath = findfile(in{a,2}, [GL_DATPATH, in{a,1}], '.abf');
    [dat, si, h] = abfload(fpath);
    
    N = size(dat, 1);
    tt = (0:N-1) .* (si * 1e-6);
    
    % plot the data
    window = [in{a,4}-0.030, in{a,5}+0.050];
    idx = (tt >= window(1)) & (tt < window(2));
    figure
    plot(tt(idx), permute(dat(idx,1,:), [1,3,2]))
    axis tight
    xlabel('time')
    ylabel('Current')
    title([in{a,1}, ' ', in{a,2}])
    
    % extract the average current in the analysis window
     window = [in{a,4}, in{a,5}]; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     Im_peak = mean(tmp, 1);
     
     % deal with the baseline response
     window = [1.32, 1.3275]; % seconds from sweep start
     idx = (tt >= window(1)) & (tt < window(2));
     tmp = permute(dat(idx,1,:), [1,3,2]);
     Im(:,a) =Im_peak - mean(tmp, 1);
    
     % customize the voltage relative to the holding potential used in the
     % experiments
     Vcmd(:,a) = linspace(-40,20,13) + in{a,3};
    
end


figure
plot(Vcmd, Im, 'linewidth', 2)
legend({in{:,2}}, 'location', 'northwest')
xlabel('Voltage (mv)')
ylabel('Current (pA)')







  
  