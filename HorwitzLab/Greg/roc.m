function auc = roc(noise,sig)

data = nan(length(noise), 1);
for i = 1:length(noise)
    data(i) = sum(sig>noise(i))+0.5*sum(sig==noise(i));
end
auc = sum(data)./(length(sig)*length(noise));
