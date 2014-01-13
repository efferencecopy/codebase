
% A little demo to see how my algorithm looks in 2-D
% normalization:
% 0) Fix marginal variances to equal '1'
% 1) Fix largest marginal variance to equal '1' (squeezing along the axes)
% 2) Fix largest conditional variance to equal '1' (squeezing along the principal axes)

function AdaptCovTest2()
    kernel = [1 .2];
    kernel = kernel./norm(kernel);
    ntrials = 500;
    nstimpertrial = 100;
    prewhiten = 0;
    alpha = 0.15;
    normalization = 2;
    longaxes = [];
    niters = 25;
    figure;
    for j = 1:niters
        StimCov = eye(2);
        nspikes = 0;
        subplot(ceil(sqrt(niters)),ceil(sqrt(niters)),j);
        hold on;
        plot([0 kernel(1)],[0 kernel(2)]);
        for i = 1:ntrials
            x = normrnd(0,1,2, nstimpertrial);
            if (prewhiten)
                [u,s,v] = svd(StimCov);
                wz = v*sqrt(s);
                iwz = inv(wz);
                M = wz*u';  % same as M = sqrtm(StimCov)
            else
                M = sqrtm(StimCov);
            end
            stim = M*x;
            spikes = modelneuron(stim, kernel);
            if (sum(spikes) == 0)
                continue
            end
            nspikes = nspikes+sum(spikes);
            STstim = stim(:,logical(spikes));
            outerprod = STstim*STstim';
            if (prewhiten)
                StimCov_transformed = (iwz*StimCov*iwz')+(iwz*(alpha*outerprod/sum(spikes))*iwz');
                StimCov = wz*StimCov_transformed*wz';  % back to the stimulus space
            else
                StimCov = StimCov+(alpha*outerprod/sum(spikes));
            end
            % Normalizing variances
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
                StimCov = v*inv(sqrt(d(1)))*v'*StimCov*v*inv(sqrt(d(1)))*v';
            end
            
            c = cos(linspace(-pi,pi,100));
            s = sin(linspace(-pi,pi,100));
            y = sqrtm(StimCov)*[c;s];
            col = [(i-1)/(ntrials-1) 1-(i-1)/(ntrials-1) 0];
            plot(y(1,:),y(2,:),'Color',col);
        end
        [v,d] = eig(StimCov);
        idx = find(diag(d) == max(diag(d)));
        v = v*sign(kernel*v(:,idx));
        longaxes = [longaxes, v(:,idx)];
        axis equal;
        set(gca,'XTick',[],'YTick',[]);
        drawnow;
    end

    thetas = atan2(longaxes(2,:),longaxes(1,:));
    figure; axes; hold on;
    hist(thetas,linspace(-pi,pi,30));
    plot(atan2(kernel(2),kernel(1)),0,'m*','MarkerSize',20);
    set(gca,'XLim',[-pi pi]);
   % keyboard
    
    % Support function that simulates a neuron
    function spikes = modelneuron(stim,kernel)
        projs = kernel*stim;
        p = 1./(1+exp(-projs+1));
        %p = 1./(1+exp(-abs(projs)+1));
        spikes = binornd(1, p);
    end
end