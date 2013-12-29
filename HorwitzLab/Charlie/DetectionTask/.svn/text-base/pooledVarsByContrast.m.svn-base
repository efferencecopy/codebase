function poolVars = pooledVarsByContrast(synthPopVars, rho)
        %takes a matrix of nNeurons X mContrasts and converts it into a
        %matris of 1PooledNeuron x mContrasts. Each element in
        %"synthPopVars is the variance of a single neuron at a particular
        %contrast. I'll solve for the cov mtx at each contrast and then sum
        %all the elements to get the population variance at each contrast
        %
        %  ASSUMES SUMS OF RVS !!!!!!
        %  DIVIDE THE OUTPUT BY POOLSIZE^2 FOR AVGS OF RVs
        %
        
        popSize = size(synthPopVars,1);
        if popSize <= 1;
            poolVars = synthPopVars;
            return %the cov mtx will just be the variance itself
        end
        
        poolVars = nan(1,size(synthPopVars,2));
        for cntrst = 1:size(synthPopVars,2);
            
            %assemble the cov mtx
            covMtx = nan(popSize, popSize);
            for row = 1:popSize;
                for col = 1:popSize;
                    if row == col;
                        covMtx(row,col) = synthPopVars(row, cntrst);
                    else
                        %the covariance of 2 RVs is the correlation times
                        %the product of the SDs
                        covMtx(row,col) = rho .* sqrt(synthPopVars(row,cntrst).*synthPopVars(col,cntrst));
                    end
                end
            end
            
            %sum the covMtx to obtain pooled variance for the sum of RVs
            poolVars(cntrst) = sum(covMtx(:));
        end
end