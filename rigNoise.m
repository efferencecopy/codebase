 %% Load in the data files and make a PS

 fin
 
 % some stuff for the analysis
 snpLength = 1;     % in sec
 idx = 0;
 fileName = 'withAndWithoutIsolator';
 
 % go to the data directory and pull out the available file names
 cd(GL_DATPATH);
 cd('../')
 cd(['./Rig_Noise', filesep, fileName])
 d = dir;

 
 for a = 1:numel(d)
     if isempty(strfind(d(a).name, '.abf'))
         continue
     end
     
     %open the data file
     fprintf('file %d of %d\n', a, numel(d));
     [dat, si, h] = abfload(d(a).name);
     
     % partition the data into snips
     tmp = dat(:,2);
     sampRate = 1./(si.*10.^(-6));
     nSamps = floor(snpLength .* sampRate);
     nSnips = floor(numel(tmp)/nSamps);
     tmp = reshape(tmp(1:(nSnips.*nSamps)), nSamps, nSnips); % time goes down the columns.
     
     
     % make the freq axis
     N = nSamps;
     if rem(N,2)
         k = -((N-1)./2):((N-1)./2); % pos freqs when N is Odd. For two sided freqs: k = -((N-1)./2):((N-1)./2)
     else
         k = -(N./2):((N./2)-1); % 0:(N./2) => abs(neg Freqs) when N is even. For two sided freqs: k = -(N./2):((N./2)-1)
     end
     freqAxis = (k./N).*sampRate;
     l_pos = freqAxis>0;
     
     
     
     % compute the average PS
     FFT_coeffs = fft(tmp, [], 1);
     FFT_mag = abs(FFT_coeffs)./N;
     PS = (FFT_mag.^2) .* 2;  % technically I need to ignore neg freq b/c i multiplied by 2
     PS = mean(PS, 2);
     PS = fftshift(PS);
     dF = freqAxis(2) - freqAxis(1);
     PowPerHZ = PS(l_pos) .* dF;
     
     % store the PS
     idx = idx+1;
     PS_avg{idx} = PowPerHZ;
     ff{idx} = freqAxis(l_pos);
     leg{idx} = d(a).name;
    
 end
 
 
 
 
 %Plot the PSs
 figure; cmap = colormap('jet'); close;
 inds = round(linspace(1, size(cmap,1), numel(PS_avg)));
 clrs = cmap(inds,:);

 
 figure, hold on,
 set(gcf, 'position', [77           5        1221         801])
 for a = 1:numel(PS_avg)
     plot(ff{a}, PS_avg{a}, 'color', clrs(a,:));
     set(gca, 'xscale', 'log', 'yscale', 'log')
 end
 axis tight
 legend(leg)
 ylabel('pA^2 per Hz')
 
 % PS individually
 figure
 for a = 1:numel(PS_avg)
     subplot(1, numel(PS_avg), a)
     loglog(ff{a}, PS_avg{a})
     axis tight
     title(leg{a})
 end
 set(gcf, 'position', [ -472         584        3963         205])
 
ymax = 0;
ymin = inf;
 for a = 1:numel(PS_avg)
     subplot(1,numel(PS_avg), a)
     yvals = get(gca, 'ylim');
     if yvals(1) < ymin; ymin = yvals(1); end
     if yvals(2) > ymax; ymax = yvals(2); end
 end
 
 for a = 1:numel(PS_avg)
     subplot(1, numel(PS_avg), a)
     ylim([ymin, ymax])
     xlim([1 ff{a}(end)])
 end
 
 
 