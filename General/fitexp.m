function fit = fitexp(in, bkgnd)

%
% FUNCTION fitexp
%
%  example: fit_fxnHand = fitexp(rawData)
%
% Should accept rawdata inputs. Fits the sum of three exponentials of the
% form:
%
% fit = [A_1 * exp(-tt/tau_1)] + [A_2 * exp(-tt/tau_2)] + [A_3 * exp(-tt/tau_3)] + C
%
% Returns a function handle that accepts a tt (time) vector, which will
% create the smooth-fitted curve. Times are with refrence to the beginning
% of the raw data set (such that t=0 is the first index to "rawData"
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
tt = [0:(numel(raw)-1)]';

% do the fit
options.MaxIter = 100000;
options.MaxFunEvals = 100000;
guesses = [raw(1).*0.2, abs(raw(1)), 10,...
           raw(1).*0.5, abs(raw(1)), 10,...
           raw(1)     , abs(raw(1)), 10, bkgnd]; 
params = fminsearch(@exp_sse, guesses, options);



% plot the result
fit = (params(1) * exp((-tt+params(3))/params(2))) + ...
      (params(4) * exp((-tt+params(6))/params(5))) + ...
      (params(7) * exp((-tt+params(9))/params(8))) + ...
       params(10);
          
          
figure
subplot(2,1,1)
hold on,
plot(raw, 'b', 'linewidth', 2)
plot(fit, 'r', 'linewidth', 2)
hold off
subplot(2,1,2)
hold on,
plot(raw-fit)
plot([0 numel(raw)], [0 0], 'k', 'linewidth', 2)
set(gcf, 'position', [846    19   596   787]);
drawnow

    function sse = exp_sse(input)
        fit = (input(1) * exp((-tt+input(3))/input(2))) + ...
              (input(4) * exp((-tt+input(6))/input(5))) + ...
              (input(7) * exp((-tt+input(9))/input(8))) + ...
              input(10);

        sse = sum((raw - fit).^2);
        
% %         % figure out the asympotic value. If bigger than the bkgnd value,
% %         % than set the SSE to INF
% %         tmp = 16e3; % twice the IPI for a 5 Hz stim, assuming 40kHz sampRate
% %         fit_asmy = (input(1) * exp((-tmp+input(3))/input(2))) + ...
% %                    (input(4) * exp((-tmp+input(6))/input(5))) + ...
% %                    (input(7) * exp((-tmp+input(9))/input(8))) + ...
% %                     input(10);
% %                 
% % %         if fit_asmy > bkgnd
% % %             sse = Inf;
% % %         end
        if any(input([1,4,7])>=0) % for EPSCs make sure the exponential starts negative
            sse = Inf;
        end
        if any(input([2,5,8])<=0) % the time constants must be positive.
            sse = Inf;
        end


       
    end

end


