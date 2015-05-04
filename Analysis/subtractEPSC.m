function ipsc = subtractEPSC(epsc, mpsc, method)


% baseline subtract all the traces
base_epsc = mean(epsc.Im(epsc.tvec<0));
base_mpsc = mean(mpsc.Im(mpsc.tvec<0));

epsc.Im = epsc.Im - base_epsc;
mpsc.Im = mpsc.Im - base_mpsc;

%sanity check
figure, hold on,
plot(epsc.tvec, epsc.Im, 'b')
plot(mpsc.tvec, mpsc.Im, 'r')


% down sample
N = 5;
downSamp_mpsc = decimate(mpsc.Im, N);
downSamp_tvec = decimate(mpsc.tvec, N);



method = 'peak'; % 'peak', or 'didt'
switch method
    case 'didt'
        
        % look at the derivitive of the Im for the down-sampled psc. Down
        % sample b/c di goes to zero when the kinetics are slow relative to the
        % sampling rate. Find where the di/dt is maximally negative following
        % the first pulse. This is where the mpsc transitions from linear to
        % compressive.
        didt = [0,diff(downSamp_mpsc)];
        [~, first_pulse_idx] = min(abs(downSamp_tvec - mpsc.tcross(1)));
        [~, second_pulse_idx] = min(abs(downSamp_tvec - mpsc.tcross(2)));
        [maxSlope_val, maxSlope_idx] = min(didt(first_pulse_idx:second_pulse_idx));
        maxSlope_idx = maxSlope_idx + first_pulse_idx- 1;
        assert(maxSlope_val == didt(maxSlope_idx), 'ERROR: max slope value is incorrect')
        maxSlope_time = downSamp_tvec(maxSlope_idx);
        
        % more sanity check
        figure, hold on,
        plot(downSamp_tvec, downSamp_mpsc, 'r')
        plot(downSamp_tvec, didt, 'g')
        plot(downSamp_tvec(first_pulse_idx), didt(first_pulse_idx), 'om', 'markersize', 10, 'markerfacecolor', 'm')
        plot(downSamp_tvec(second_pulse_idx), didt(second_pulse_idx), 'om', 'markersize', 10, 'markerfacecolor', 'm')
        plot(downSamp_tvec(maxSlope_idx), didt(maxSlope_idx), 'oc', 'markersize', 10, 'markerfacecolor', 'c')
        
        % translate back to the normal sampling rate.
        [~, window_end_idx] = min(abs(mpsc.tvec - maxSlope_time));
        [~, window_start_idx] = min(abs(mpsc.tvec - mpsc.tcross(1)));
        
        % pull out the snippets to scale
        snippet_epsc = epsc.Im(window_start_idx:window_end_idx);
        snippet_mpsc = mpsc.Im(window_start_idx:window_end_idx);
        
    case 'peak'
        
        % find the peak in the down sampled trace
        [~, first_pulse_idx] = min(abs(downSamp_tvec - mpsc.tcross(1)));
        [~, second_pulse_idx] = min(abs(downSamp_tvec - mpsc.tcross(2)));
        peak_val = min(downSamp_mpsc(first_pulse_idx:second_pulse_idx));
        Im_start = 0.1 .* peak_val;
        Im_end = 0.85 .* peak_val;
        
        % now go back to the upsampled time world, and figure out when the
        % mpsc crosses 5% and 85% of the peak value
        peak_idx = find(downSamp_mpsc == peak_val, 1);
        peak_time = downSamp_tvec(peak_idx);
        t_window = (mpsc.tvec>mpsc.tcross(1)) & (mpsc.tvec<peak_time);
        valid = (mpsc.Im <= Im_start) & (mpsc.Im >= Im_end) & t_window;
        
        % plot the rise time + the "valid" points
        figure, hold on,
        plot(mpsc.tvec(t_window), mpsc.Im(t_window), 'r')
        plot(mpsc.tvec(valid), mpsc.Im(valid), 'k.')
        
        % pull out the snippets to scale
        snippet_epsc = epsc.Im(valid);
        snippet_mpsc = mpsc.Im(valid);
end


% find the scale factor
ep = snippet_epsc;
mp = snippet_mpsc;
minFxn = @(c) sum((ep.*c - mp).^2);
scale = fminsearch(@(x) minFxn(x), 1);

% plot the original snippets, and the scaled epsc
figure, hold on,
plot(snippet_epsc, 'b')
plot(snippet_mpsc, 'r')
plot(snippet_epsc.*scale, 'c.')

% subtract off the scalled epsc, and plot
ipsc = mpsc.Im - epsc.Im.*scale;


% plot everything together
figure, hold on,
plot(mpsc.tvec, ipsc, 'r');
plot(epsc.tvec, epsc.Im.*scale, 'b');
plot(mpsc.tvec, mpsc.Im, 'm')





