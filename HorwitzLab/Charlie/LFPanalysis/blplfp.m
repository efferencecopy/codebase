function blp = blplfp(stro, lfpidx, band)

    % pull out the raw lfps
    sampFreq = stro.sum.analog.storeRates{1};
    rawLFP = stro.ras(:, lfpidx);
    
    
    nTrials = length(rawLFP);
   
    %setup a notch filter:
    flow = 59;
    fhigh = 61;
    order = 4;
    r =.01;
    nyquist = sampFreq/2;
    Wn = [(flow/nyquist) , fhigh/nyquist];
    [B, A] = cheby1(order/2,r, Wn, 'stop');
    Bcell = mat2cell(repmat(B,nTrials,1), ones(nTrials,1));
    Acell = mat2cell(repmat(A,nTrials,1), ones(nTrials,1));
    
    if numel(band) == 2; %assume cheby filtering
        
        %run the notch filter after windowing
        windowCell = cellfun(@(x) blackman(numel(x)), rawLFP, 'uniformoutput', 0);
        windowedLFP = cellfun(@(x,y) x.*y, rawLFP, windowCell, 'uniformoutput', 0);
        windowedLFP = cellfun(@filtfilt, Bcell, Acell, rawLFP, 'uniformoutput', 0);
        
        %setup the bandpass filter
        flow = band(1);
        fhigh = band(2);
        order = 4;
        r =.01;
        Wn = [flow/nyquist , fhigh/nyquist];
        [B, A] = cheby1(order/2, r, Wn);
        
        %filter the lfps
        Bcell = mat2cell(repmat(B,nTrials,1), ones(nTrials,1));
        Acell = mat2cell(repmat(A,nTrials,1), ones(nTrials,1));
        bandLimitLFPs = cellfun(@filtfilt, Bcell, Acell, windowedLFP, 'uniformoutput', 0);
        blp = bandLimitLFPs;
        
    elseif numel(band) == 1; %assume greg's thing
        
        % run the notch filter but do not window
        notchLFP = cellfun(@filtfilt, Bcell, Acell, rawLFP, 'uniformoutput', 0);
        
        centerfreq = band; % Hz
        cyclespersample = centerfreq./sampFreq;
        ncyclesper6sigma = 6;  % Controls the bandwidth
        nsamplesper6sigma = ceil(ncyclesper6sigma/cyclespersample);
        ncycles = nsamplesper6sigma*cyclespersample;
        filtkernel1 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*cos(linspace(0,2*pi*ncycles,nsamplesper6sigma));
        filtkernel2 = normpdf(linspace(-3,3,nsamplesper6sigma),0,1).*sin(linspace(0,2*pi*ncycles,nsamplesper6sigma));
        filtkernel1 = filtkernel1./norm(filtkernel1);
        filtkernel2 = filtkernel2./norm(filtkernel2);
        
        %perform the convolution
        kernel1 = mat2cell(repmat(filtkernel1,nTrials,1) , ones(nTrials,1));
        kernel2 = mat2cell(repmat(filtkernel2,nTrials,1) , ones(nTrials,1));
        %blp = cellfun(@(x, k1, k2) conv(x, k1,'same').^2+conv(x,k2,'same').^2, notchLFP, kernel1, kernel2, 'uniformoutput', 0);
        blp = cellfun(@(x, k1, k2) ((conv(x, k1,'same').^2) ./ numel(x)), notchLFP, kernel1, 'uniformoutput', 0);
    end

    
    %determine which trials saturated the AD board, and thus should be
    %ignored. The oob trials are included, but the values of the anlg
    %channel are sent to NaNs.
    oob = stro.checkTrials(lfpidx);
    blp(oob) = cellfun(@(x) x.*nan, blp(oob), 'uniformoutput', 0);
end

