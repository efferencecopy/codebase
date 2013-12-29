% This is some code to test out ideas for adaptively manipulating the spectrum 
% of reverse correlation stimulus to focus our energy in the relevant subspace 
% normalization:
% 0) Fix marginal variances to equal '1'
% 1) Fix largest marginal variance to equal '1' (squeezing along the axes)
% 2) Fix largest conditional variance to equal '1' (squeezing along the principal axes)


function AdaptCovTest
    prewhiten = 0;
    normalization = 0;
    nfilters = 2;
    angle = pi/2;
    ndims = 200;
    nstimpertrial = 300;
    ntrials = 300;
    StimCov = eye(ndims);
    alpha = 0.01;
    
    kernels = mkkernels(nfilters,[1:nfilters] ,angle);
  %  corrcoef(kernels');
  %  rank(kernels);
    nspikes = 0;
    hwait = waitbar(0,'Simulating experiment...');
    for i = 1:ntrials
        waitbar(i/ntrials,hwait);
        x = normrnd(0,1,ndims, nstimpertrial);
        if (prewhiten)
            [u,s,v] = svd(StimCov);
            wz = v*sqrt(s);
            iwz = inv(wz);
            M = wz*u';  % same as M = sqrtm(StimCov)
        else
            M = sqrtm(StimCov);
        end
        stim = M*x;
        spikes = modelneuron(stim, kernels);
        nspikes = nspikes+sum(spikes);
        if (sum(spikes) == 0)
            continue
        end
        
        STstim = stim(:,logical(spikes));
        outerprod = STstim*STstim';
        if (prewhiten)
            StimCov_transformed = (iwz*StimCov*iwz')+(iwz*(alpha*outerprod/sum(spikes))*iwz');
            StimCov = wz*StimCov_transformed*wz';  % back to the stimulus space
        else
            StimCov = StimCov+(alpha*outerprod/sum(spikes));
        end
        if (normalization == 0)
            sds = sqrt(diag(StimCov));
            normmat = (1./sds)*(1./sds)';
            StimCov = StimCov.*normmat;
        elseif(normalization == 1)
            sds = sqrt(diag(StimCov));
            sds = max(sds);
            normmat = (1./sds)*(1./sds)';
            StimCov = StimCov.*normmat;
        else %  normalization == 3

            [u,d,v] = svd(StimCov);
            
            %d = d(size(StimCov,1)-nspikes+1:end,size(StimCov,1)-nspikes+1:end);
            %v = v(:,size(StimCov,1)-nspikes+1:end);
            %normmat = v*diag(1./sqrt(diag(d)))*v';
            %normmat = inv(v*sqrt(d));
            greg = v*inv(sqrt(d(1)))*v'*StimCov*v*inv(sqrt(d(1)))*v';
            StimCov = greg;
        end 
    end
    close(hwait);
    [v,d] = eig(StimCov);
    [sortd,idx] = sort(diag(d));

    figure;
    subplot(2,3,1);
    plot(diag(StimCov),'k.')  % See if stimulus variance is growing or shrinking
    subplot(2,3,2);
    plot(sortd,'k.');
    subplot(2,3,3);
    plot(v(:,idx(end)));
    subplot(2,3,4);
    hold on;
    c = ['r','g','b','m','k'];
    for j = [1:nfilters]
        plot(v(:,idx(end-j+1)),'Color',c(j));
    end
    subplot(2,3,5);
    plot(kernels');
    %imagesc(StimCov);
    subplot(2,3,6);
    set(gca,'Visible','off')
    text(0,1,['nspikes = ',num2str(nspikes)]);
    text(0,.8,['condn number = ',num2str(cond(StimCov))]);
 
    
    % Support function that creates kernels
    function kernels = mkkernels(nfilters,freq,angle)
        envelope = normpdf(linspace(-5,5,ndims),0,1);
        for i = 1:nfilters
            kernel = envelope.*cos(freq(i).*linspace(-5,5,ndims)+(i-1)*angle);
            kernels(i,:) = kernel./norm(kernel);
        end
    end

    % Support function that simulates a neuron
    function spikes = modelneuron(stim,kernels)    
        nfilters = size(kernels,1);
        projs = [];
        for counter = 1:nfilters
            projs(counter,:) = kernels(counter,:)*stim;
        end
        if (nfilters == 1)
            p = 1./(1+exp(-projs+1));
            p = 1./(1+exp(-abs(projs)+1));
        elseif (nfilters == 2)
            p(1,:) = 1./(1+exp(-projs(1,:)+1));
            p(2,:) = 1./(1+exp(-abs(projs(2,:))+1));
            p = prod(p);
        end
        if (nfilters > 2)
            p = 1./(1+exp(-abs(projs)+1));
            p = sum(p);
            p = p./max(p);
        end
        spikes = binornd(1, p);
        if any(isnan(spikes))
            keyboard
        end
    end
end

    