function hand = fillinterval(x,mean,bottom,top,color)
% FILLINTERVAL plots interval estimate of the mean. The function only plots
% the interval given the input data, it does not compute it.
%
% The inverval is displayed as a filled area.

% make everything a column
x      = x(:);
mean   = mean(:);
bottom = bottom(:);
top    = top(:);

% plot the filled area
f = fill([x; flip(x)],...
     [bottom; flip(top)],...
     color);
set(f,'FaceAlpha',0.1)
set(f,'EdgeAlpha',0.1)

% add the line (mean)
hold on
if nargout > 0
    hand = plot(x,mean,'color',color,'linewidth',1);
else
    plot(x,mean,'color',color,'linewidth',1)
end

end