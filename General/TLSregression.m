function [slope, intercept] = TLSregression(X,Y)


% preform total least squares regression via PCA
eigvecs = pca([X,Y]);
PC1 = eigvecs(:,1);

% solve for 2D slope and intercept
slope = PC1(2)./PC1(1);
pointVal = mean([X,Y],1);
intercept = pointVal(2) - (slope.*pointVal(1));
