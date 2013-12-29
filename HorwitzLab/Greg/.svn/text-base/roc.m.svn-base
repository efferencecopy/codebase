function out = roc(x,y)
data = [];
for i = 1:length(x)
    data(i) = sum(y>x(i))+0.5*sum(y==x(i));
end
out = sum(data)./(length(y)*length(x));
end