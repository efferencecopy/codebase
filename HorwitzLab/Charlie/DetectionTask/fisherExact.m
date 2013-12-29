function [h, p] = fisherExact(f11,f21,f12,f22,alpha)

%  FISHER EXACT TEST
%
%    EXAMPLE [h, p] = fisherExact(f11,f21,f12,f22,alpha);
%
% Calcultes the probablity that the marginals of a 2x2 contingency tables
% are independent distributions. The input arguments are the counts data
% that are entered into the cells of the table (i.e., f11 is the count
% associated with the first row and column). If alpha is not specified, the
% default alpha = 0.05; 
%
% CAH 2/11
%
% See Biostatistical Analysis (4th Ed). Jarrol Zar. pp 543-552

if ~exist('alpha', 'var')
    alpha = 0.05;
end

N = sum([f11,f21,f12,f22]);
R1 = f11 + f12;
R2 = f21 + f22;
C1 = f11 + f21;
C2 = f12 + f22;
numerator = log10(factorial(R1)) + log10(factorial(R2)) + log10(factorial(C1)) + log10(factorial(C2)) - log10(factorial(N));
denominator = log10(factorial(f11)) + log10(factorial(f12)) + log10(factorial(f21)) + log10(factorial(f22));
pCrit =  10^(numerator - denominator);


%now find the probabilities associated with all possible table outcomes
%(whos marginal totals are the same as the emperical table's).
phat = nan(C1+1,1);
for i = 0:C1;
    
    f11hat = i;
    f12hat = R1 - f11hat;
    f21hat = C1 - f11hat;
    f22hat = C2 - f12hat;
    
    if any([f11hat, f12hat, f21hat, f22hat] < 0)
        continue %i.e., some of the table entries are negative...
    end
    
    denominator = log10(factorial(f11hat)) + log10(factorial(f12hat)) + log10(factorial(f21hat)) + log10(factorial(f22hat));
    phat(i+1) = 10^(numerator - denominator);
end

p = sum(phat(phat<=pCrit));
h = p < alpha;