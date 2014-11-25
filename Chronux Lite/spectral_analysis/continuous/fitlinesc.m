function [datafit,Amps,freqs,Fval,sig]=fitlinesc(data,params,p,plt,f0)
% fits significant sine waves to data (continuous data).
%
% Usage: [datafit,Amps,freqs,Fval,sig]=fitlinesc(data,params,p,plt,f0)
%
%  Inputs:  
% Note that units of Fs, fpass have to be consistent.
%       data        (data in [N,C] i.e. time x channels/trials or a single
%       vector) - required.
%       params      structure containing parameters - params has the
%       following fields: tapers, Fs, fpass, pad
%           tapers : precalculated tapers from dpss or in the one of the following
%                    forms: 
%                   (1) A numeric vector [TW K] where TW is the
%                       time-bandwidth product and K is the number of
%                       tapers to be used (less than or equal to
%                       2TW-1). 
%                   (2) A numeric vector [W T p] where W is the
%                       bandwidth, T is the duration of the data and p 
%                       is an integer such that 2TW-p tapers are used. In
%                       this form there is no default i.e. to specify
%                       the bandwidth, you have to specify T and p as
%                       well. Note that the units of W and T have to be
%                       consistent: if W is in Hz, T must be in seconds
%                       and vice versa. Note that these units must also
%                       be consistent with the units of params.Fs: W can
%                       be in Hz if and only if params.Fs is in Hz.
%                       The default is to use form 1 with TW=3 and K=5
%
%	        Fs 	        (sampling frequency) -- optional. Defaults to 1.
%               fpass       (frequency band to be used in the calculation in the form
%                                   [fmin fmax])- optional. 
%                                   Default all frequencies between 0 and Fs/2
%	        pad		    (padding factor for the FFT) - optional (can take values -1,0,1,2...). 
%                    -1 corresponds to no padding, 0 corresponds to padding
%                    to the next highest power of 2 etc.
%			      	 e.g. For N = 500, if PAD = -1, we do not pad; if PAD = 0, we pad the FFT
%			      	 to 512 points, if pad=1, we pad to 1024 points etc.
%			      	 Defaults to 0.
%	    p		    (P-value to calculate error bars for) - optional. 
%                           Defaults to 0.05/N where N is data length.
%       plt         (y/n for plot and no plot respectively) - plots the
%       Fratio at all frequencies if y
%       f0          frequencies at which you want to remove the
%                   lines - if unspecified the program
%                   will compute the significant lines
%
%
%  Outputs: 
%       datafit        (linear superposition of fitted sine waves)
%       Amps           (amplitudes at significant frequencies)
%       freqs          (significant frequencies)
%       Fval           (Fstatistic at all frequencies)
%       sig            (significance level for F distribution p value of p)


data=change_row_to_column(data);
[N,C]=size(data);

if nargin < 2 || isempty(params); params=[]; end;

[tapers,~,Fs,~,~,~,params]=getparams(params);
clear pad fpass err trialave;

if nargin < 3 || isempty(p);p=0.05/N;end;
if nargin < 4 || isempty(plt); plt='n'; end;
if nargin < 5; f0=[]; end;

params.tapers=dpsschk(tapers,N,Fs); % calculate the tapers
[Fval,A,f,sig] = ftestc(data,params,p,plt);

if isempty(f0);
   
   [~, loc] = findpeaks(Fval,'threshold', sig);
   
   freqs=f(loc);
   
   Amps=A(loc);
   datafit=exp(1i*2*pi*(0:N-1)'*freqs/Fs)*A(loc)+exp(-1i*2*pi*(0:N-1)'*freqs/Fs)*conj(A(loc));
   
else
    indx = zeros(length(f0), 1);
    for n=1:length(f0);
        [~,indx(n)]=min(abs(f-f0(n)));
    end;
    freqs=f(indx);
    Amps = A(indx);
    datafit=exp(1i*2*pi*(0:N-1)'*freqs/Fs)*A(indx)+exp(-1i*2*pi*(0:N-1)'*freqs/Fs)*conj(A(indx));

end




