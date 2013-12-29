function f = HalfSquareFit(b, x, response)

% function f = HalfSquareFit(b, x, response)
% Fits a half-squaring contrast response function (with threshold)
% to data asuming Poisson error. 
%
% One of two different models is fit, depending on the length of the
% parameter vector, b. If it has 3 elements, this model is fit:
%
% prediction = b(1)+b(3)*(max(x-b(2),0).^2);
%
% and 'x' is assumed to be a vector of contrasts
%
% If b has four elements, x is assumed to have two columns. The first 
% is the vector of contrasts, the second is a binary vector of 'indicator'
% variables (I). This model is fit:
%
%       for I = 0: pred = b(1)+b(3)*(max(x-b(2),0).^2);
%       for I = 1: pred = b(1)+b(3)*(max(x.*b(4)-b(2),0).^2);
try
if (length(b) == 3)
    fittype = 'simple';
    % b(1)+b(3)*max(x-b(2),0)^2
else
    fittype = 'yoked';
    L = logical(x(:,2));
end

if (size(x,2) > 1 && strcmp(fittype,'simple'))
    error('too many predictors for simple fit');
end
if (size(x,2) ~=2 && ~strcmp(fittype,'simple'))
    error('expecting 2 columns for design matrix');
end

if (any(isnan(b)))
    f = Inf;
else
    if (strcmp(fittype,'simple'))
        prediction = (b(1)+b(3)*(max(x-b(2),0).^2))+eps;
        f = sum(prediction)-(log(prediction')*response);  % -1 * log-likelihood (Poisson)
    end
    if (strcmp(fittype,'yoked'))
        pred1 = (b(1)+b(3)*(max(x(~L)-b(2),0).^2))+eps;
        pred2 = (b(1)+b(3)*(max(x(L).*b(4)-b(2),0).^2))+eps;
        f1 = sum(pred1)-(log(pred1')*response(~L));
        f2 = sum(pred2)-(log(pred2')*response(L));
        f = f1+f2;
    end
end
catch
    disp('in HalfSquareFit_old')
    keyboard
end
    
    