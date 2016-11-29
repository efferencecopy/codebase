function [top_int, bot_int, X_mod, Y_mod] = regression_line_ci(alpha,beta,x,y,varargin)


%[TOP_INT, BOT_INT] = REGRESSION_LINE_CI(ALPHA,BETA,X,Y)
%creates two curves marking the 1 - ALPHA confidence interval for the
%regression line, given BETA coefficience (BETA(1) = intercept, BETA(2) =
%slope). This is the format of STATS.beta when using 
%STATS = REGSTATS(Y,X,'linear','beta');
%[TOP_INT, BOT_INT] = REGRESSION_LINE_CI(ALPHA,BETA,X,Y,N_PTS) defines the
%number of points at which the funnel plot is defined. Default = 100
%[TOP_INT, BOT_INT] = REGRESSION_LINE_CI(ALPHA,BETA,X,Y,N_PTS,XMIN,XMAX)
%defines the range of x values over which the funnel plot is defined

if(length(x) ~= length(y))
    error(message('regression_line_ci:x and y size mismatch')); 
end

% strip out nans from either X or Y
l_nan = isnan(x) | isnan(y);
x = x(~l_nan);
y = y(~l_nan);

N = length(x);


x_min = min(x);
x_max = max(x);
n_pts = 100;

if(nargin > 4)
    n_pts = varargin{1};
end
if(nargin > 6)
    x_min = varargin{2};
    x_max = varargin{3};
end

X_mod = x_min:(x_max-x_min)/n_pts:x_max;
Y_mod = ones(size(X_mod))*beta(1) + beta(2)*X_mod;

SE_y_cond_x = sum((y - beta(1)*ones(size(y))-beta(2)*x).^2)/(N-2);
SSX = (N-1)*var(x);
SE_Y = SE_y_cond_x*(ones(size(X_mod))*(1/N + (mean(x)^2)/SSX) + (X_mod.^2 - 2*mean(x)*X_mod)/SSX);

Yoff = (2*finv(1-alpha,2,N-2)*SE_Y).^0.5;


% SE_b0 = SE_y_cond_x*sum(x.^2)/(N*SSX)
% sqrt(SE_b0)

top_int = Y_mod + Yoff;
bot_int = Y_mod - Yoff;

