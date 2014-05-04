function [fit, startIdx] = fitexp(in, bkgnd, params)

%
% FUNCTION fitexp
%
%  example: fit_fxnHand = fitexp(rawData)
%
% Should accept rawdata inputs. Fits the sum of three exponentials of the
% form:
%
% fit = [A_1 * exp(-xx/tau_1)] + [A_2 * exp(-xx/tau_2)] + [A_3 * exp(-xx/tau_3)] + C
%
% Returns a function handle that accepts an xx (sample number) vector, which will
% create the smooth-fitted curve. Time samples are with refrence to the beginning
% of the raw data set (such that x=0 is the first index to "rawData"
%
%
% C.Hass


% Need to figure out if I should penalize for asympotic value
%
% Need to figure out if Inf is better than sse*10 for penalty
%
% Need to figure out how to increase fit quality for very small PSCs
%
% Need to adapt the function for EPSCs and IPSCs.


% Try to find a reasonable starting place by looking at the 1nd derivitive.
% When this is maximally positive, than for a EPSC, this is the time of
% maximum positive slope, and is when the decaying exponential take over
% from the other temporal dynamics around the peak inward current.
dxdt = diff(in);
windowSize = 40; % in points
dxdt = filtfilt(ones(1,windowSize)/windowSize, 1, dxdt);
[~, startIdx] = max(dxdt(windowSize:end-windowSize));
startIdx = startIdx+windowSize;

% define the snippet to fit based off the startIdx
raw = in(startIdx:end);
xx = [0:(numel(raw)-1)]';


% Make different guesses if the PSC is very small vs. the case
% where the PSC is large.
options.MaxIter = 100000;
options.MaxFunEvals = 100000;
if mean(raw(1:10)) < -20;
    guesses = [raw(1).*0.2, abs(raw(1)), 10,...
               raw(1).*0.5, abs(raw(1)), 10,...
               raw(1)     , abs(raw(1)), 10];
else
    guesses = [raw(1).*0.01, abs(raw(1)), 10,...
               raw(1).*0.2, abs(raw(1)), 10,...
               raw(1)     , abs(raw(1)), 10];
end

% do the fit   
[coeffs, ~, exitflag] = fminsearch(@exp_sse, guesses, options);



% construct the return argument
fit = @(yy) (coeffs(1) * exp((-yy+coeffs(3))/coeffs(2))) + ...
      (coeffs(4) * exp((-yy+coeffs(6))/coeffs(5))) + ...
      (coeffs(7) * exp((-yy+coeffs(9))/coeffs(8)));
          

  if isfield(params, 'debug') && params.debug
      figure
      subplot(2,1,1)
      hold on,
      plot(raw, 'b', 'linewidth', 2)
      plot(fit(xx), 'r', 'linewidth', 2)
      hold off
      title(num2str(coeffs))
      subplot(2,1,2)
      hold on,
      plot(raw-fit(xx))
      plot([0 numel(raw)], [0 0], 'k', 'linewidth', 2)
      set(gcf, 'position', [689    19   685   787]);
      drawnow
  end
  


%
% NESTED SUBFUNCTION
%
%%%%%%%%%%%%%%%%
    function sse = exp_sse(input)
        fit = (input(1) * exp((-xx+input(3))/input(2))) + ...
              (input(4) * exp((-xx+input(6))/input(5))) + ...
              (input(7) * exp((-xx+input(9))/input(8)));

        sse = sum((raw - fit).^2);
        
        % figure out the asympotic value. If bigger than the bkgnd value,
        % than set the SSE to INF
        tmp = 16e3; % twice the IPI for a 5 Hz stim, assuming 40kHz sampRate
        fit_asmy = (input(1) * exp((-tmp+input(3))/input(2))) + ...
                   (input(4) * exp((-tmp+input(6))/input(5))) + ...
                   (input(7) * exp((-tmp+input(9))/input(8)));
                
        if fit_asmy > bkgnd
            sse = Inf;
        end
        if any(input([1,4,7])>=0) % for EPSCs make sure the exponential starts negative
            sse = Inf;
        end
        if any(input([2,5,8])<=0) % the time constants must be positive.
            sse = Inf;
        end


       
    end


end


