clear
close all

% load data
load('inpainting_demo_01.mat')

% format of the first graph
stacked = true;
%#ok<*UNRCH>

%% time-domain solutions
figure
for method = 1:length(methods)
    plot(results.signals{method}(:,end))
    hold on
end
plot(y)
xlim([1 length(y)])
legend([methodnames,'original'])
set(gcf,'windowstate','maximized')

%% objective functions
figure
tl = tiledlayout(2,3);
data = NaN(maxit,length(methods));
for method = 1:length(methods)
    data(1:length(results.objectives{method}),method) = results.objectives{method};
end
if stacked
    s = tiledlayout(tl,length(methods),1,'tilespacing','tight','padding','tight');
    s.Layout.TileSpan = [2 1];
    C = colororder;
    for method = 1:length(methods)
        nexttile(s)
        plot(data(:,method),'color',C(method,:))
        legend(methodnames{method})
        ylim('padded')
        if method == 3
            ylabel('-log likelihood')
        end
    end
else
    nexttile(tl,[2 1])
    plot(data)
    ylabel('-log likelihood (may differ between methods)')
    ylim('padded')
    legend(methodnames)
end
xlabel('iteration')

%% SNRs
nexttile(tl,[2 1])
for method = 1:length(methods)
    plot(results.SNRs{method})
    hold on
end
xlabel('iteration')
ylabel('SNR (dB)')
ylim('padded')
legend(methodnames,'location','southeast')

%% relative change of objective
nexttile(tl)
for method = 1:length(methods)
    data = Inf(length(results.objectives{method}),1);
    for i = 2:length(results.objectives{method})
        data(i) = abs(results.objectives{method}(i)-results.objectives{method}(i-1))...
            /abs(results.objectives{method}(i-1));
    end
    semilogy(data)
    hold on
end
hold off
xlabel('iteration')
ylabel('relative change of the objective')
ylim('padded')
legend(methodnames)

%% relative norms
nexttile(tl)
for method = 1:length(methods)
    semilogy(results.relnorms{method})
    hold on
end
hold off
xlabel('iteration')
ylabel('relative change of the solution')
ylim('padded')
legend(methodnames)

set(gcf,'windowstate','maximized')