function hand = plotinterval(x,mean,bottom,top,color)
% PLOTINTERVAL plots interval estimate of the mean. The function only plots
% the interval given the input data, it does not compute it.
%
% The inverval is displayed by vertical lines.

% plot the line (mean)
if nargout > 0
    hand = plot(x,mean,'-o','color',color);
else
    plot(x,mean,'-o','color',color)
end

% add the interval for each value of x
for i = 1:length(x)
    line([x(i) x(i)],[bottom(i) top(i)],'color',0.75 + 0.25*color,'linewidth',1,'marker','_')
end

end