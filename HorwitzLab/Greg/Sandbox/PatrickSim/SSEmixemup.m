function [SSEmixmatrix] = SSEmixemup(scaled,Loog,xformmat,planeparams,nbootiter,SSEmixmatrix,q)

xyz = scaled*xformmat;

[th,ph,r] = cart2sph(xyz(~Loog,1),xyz(~Loog,2),xyz(~Loog,3));
planer = -1./(planeparams(1).*cos(ph).*cos(th)+planeparams(2).*cos(ph).*sin(th)+planeparams(3).*sin(ph));
planer = abs(planer);  % this ok?
resid = log(r)-log(planer);
   
nulldist = nan*ones(nbootiter,1); %#ok<NASGU>
wait_h = waitbar(0,'Bootstrapping...');
options = optimset('MaxFunEvals',50000,'MaxIter',50000,'TolFun',10^-6,'TolX',10^-6);
SSEs = [];
for i = 1:nbootiter
    waitbar(i/nbootiter, wait_h);
    fprintf('%u : %u \n',q,i)
    tmpresid = exp(resid(unidrnd(length(resid),[1 length(resid)])));
      
    tmpx = xyz(:,1);
    tmpy = xyz(:,2);
    tmpz = xyz(:,3);
       
    tmpx(~Loog) = planer.*tmpresid.*cos(ph).*cos(th);
    tmpy(~Loog) = planer.*tmpresid.*cos(ph).*sin(th);
    tmpz(~Loog) = planer.*tmpresid.*sin(ph);

    % Using surfacefitter2
    [~, SSE1] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),planeparams, options);
    [~, SSE2] = fminsearch(@(x) surfacefiterr2([tmpx tmpy tmpz],x, Loog),[planeparams;0;0;0], options);

    % Using NTsurfacefit
    %[planeparams, SSE1, quadparams, SSE2, xformmat] = NTsurfacefit([tmpx tmpy tmpz], Loog);

    SSEs(i,:) = [SSE1 SSE2]; %#ok<AGROW>
    [SSE1 SSE2] %#ok<NOPRT>
    SSEmixmatrix(i,q)=SSEs(i,1)./SSEs(i,2);
    
    figure(3);clf;
    set(gca,'XLim',[0 10])
    hist(SSEmixmatrix(:,q),100)
    drawnow
      
end
%title(['SSE Bootstrapping loop',num2str(q)]) %To display a histogram for every loop
title('Current SSE Bootstrapping loop')
close(wait_h);
