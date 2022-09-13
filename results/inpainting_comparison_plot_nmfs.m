clear
close all
addpath('../utils')

% plot intervals?
intervals = true;

% compute bootstraps?
bootstrap = false;

%% load data
fprintf('Loading data...\n')
load('inpainting_comparison_01.mat')

% consider only the NMF-based methods
methods = methods([1:2, end]);
methodnames = methodnames([1:2, end]);
maxit = NMF.maxit;

%% initialize the data arrays
fprintf('Processing the data...\n')
objective_abs = NaN(length(signums),length(glengths),maxit,length(methods));
objective_rel = NaN(length(signums),length(glengths),maxit,length(methods));
snr           = NaN(length(signums),length(glengths),maxit,length(methods));
pemoq         = NaN(length(signums),length(glengths),maxit,length(methods));
peaq          = NaN(length(signums),length(glengths),maxit,length(methods));
relnorm       = NaN(length(signums),length(glengths),maxit,length(methods));

%% fill the data arrays
rowcounter = 0; % counter of the signal & gap length combinations
sigcounter = 0; % counter of the signal
for signum = signums
    sigcounter = sigcounter + 1;
    gapcounter = 0;
    for glength = glengths
        
        rowcounter = rowcounter + 1;
        gapcounter = gapcounter + 1;
        
        %% read data
        for m = 1:length(methods)
            L = length(tables.(methods{m}).SNR{rowcounter});
        	snr(sigcounter,gapcounter,1:L,m)     = tables.(methods{m}).SNR{rowcounter};
            pemoq(sigcounter,gapcounter,1:L,m)   = tables.(methods{m}).PEMOQ{rowcounter};
            peaq(sigcounter,gapcounter,1:L,m)    = tables.(methods{m}).PEAQ{rowcounter}; 
            relnorm(sigcounter,gapcounter,1:L,m) = tables.(methods{m}).relnorms{rowcounter};
        end
        
    end
end

%% compute the bootstraps / stds
if intervals
    if bootstrap
        fprintf('Computing the bootstraps (this may take some time)...\n') %#ok<*UNRCH>
        addpath('../utils')
        [means_1, lowers_1, uppers_1] = bootstrap_est(snr);
        [means_2, lowers_2, uppers_2] = bootstrap_est(pemoq);
        [means_3, lowers_3, uppers_3] = bootstrap_est(peaq);
        [means_4, lowers_4, uppers_4] = bootstrap_est(log10(relnorm));
    else
        fprintf('Computing the stds...\n')
        means_1  = squeeze(mean(snr,1,'omitnan'));
        means_2  = squeeze(mean(pemoq,1,'omitnan'));
        means_3  = squeeze(mean(peaq,1,'omitnan'));
        means_4  = squeeze(mean(log10(relnorm),1,'omitnan'));
        std_1    = squeeze(std(snr,0,1,'omitnan'));
        std_2    = squeeze(std(pemoq,0,1,'omitnan'));
        std_3    = squeeze(std(peaq,0,1,'omitnan'));
        std_4    = squeeze(std(log10(relnorm),0,1,'omitnan'));
        mult     = tinv(0.975,length(signums)-1);
        lowers_1 = means_1 - mult*std_1/sqrt(length(signums));
        uppers_1 = means_1 + mult*std_1/sqrt(length(signums));
        lowers_2 = means_2 - mult*std_2/sqrt(length(signums));
        uppers_2 = means_2 + mult*std_2/sqrt(length(signums));
        lowers_3 = means_3 - mult*std_3/sqrt(length(signums));
        uppers_3 = means_3 + mult*std_3/sqrt(length(signums));
        lowers_4 = means_4 - mult*std_4/sqrt(length(signums));
        uppers_4 = means_4 + mult*std_4/sqrt(length(signums));
    end
else
    means_1 = squeeze(mean(snr,1,'omitnan'));
    means_3 = squeeze(mean(pemoq,1,'omitnan'));
    means_3 = squeeze(mean(peaq,1,'omitnan'));
    means_4 = squeeze(mean(log10(relnorm),1,'omitnan'));
end

%% plot
if intervals
    if bootstrap
        fprintf('Plotting the bootstraps...\n')
    else
        fprintf('Plotting the means + stds...\n')
    end
else
    fprintf('Plotting...\n')
end
colors = colororder;
figure('visible','off')
T  = tiledlayout(length(glengths),4,'tileindexing','rowmajor');
h  = gobjects(length(methods),1);
ax = gobjects(4,length(glengths));
for gapcounter = 1:length(glengths)
    
    % SNR
    ax(1,gapcounter) = nexttile;
    for m = 1:length(methods)
        if intervals
            h(m) = fillinterval(...
                1:maxit,...
                means_1(gapcounter,:,m),...
                lowers_1(gapcounter,:,m),...
                uppers_1(gapcounter,:,m),...
                colors(m,:));
            hold on
        else
            h(m) = plot(1:maxit,means_1(gapcounter,:,m),'color',colors(m,:));
            hold on
        end
    end
    legend(h,methodnames,'location','southeast')
    xlabel('iteration')
    ylabel('SNR (dB)')
    title(sprintf('gap length %d ms',glengths(gapcounter)))
    xlim([0 50])
    
    % PEMO-Q
    ax(2,gapcounter) = nexttile;
    for m = 1:length(methods)
        if intervals
            h(m) = fillinterval(...
                1:maxit,...
                means_2(gapcounter,:,m),...
                lowers_2(gapcounter,:,m),...
                uppers_2(gapcounter,:,m),...
                colors(m,:));
            hold on
        else
            h(m) = plot(1:maxit,means_2(gapcounter,:,m),'color',colors(m,:));
            hold on
        end
    end
    legend(h,methodnames,'location','southeast')
    xlabel('iteration')
    ylabel('PEMO-Q ODG')
    title(sprintf('gap length %d ms',glengths(gapcounter)))
    xlim([0 50])
    
    % PEAQ
    ax(3,gapcounter) = nexttile;
    for m = 1:length(methods)
        if intervals
            h(m) = fillinterval(...
                1:maxit,...
                means_3(gapcounter,:,m),...
                lowers_3(gapcounter,:,m),...
                uppers_3(gapcounter,:,m),...
                colors(m,:));
            hold on
        else
            h(m) = plot(1:maxit,means_3(gapcounter,:,m),'color',colors(m,:));
            hold on
        end
    end
    legend(h,methodnames,'location','southeast')
    xlabel('iteration')
    ylabel('PEAQ ODG')
    title(sprintf('gap length %d ms',glengths(gapcounter)))
    xlim([0 50])
    
    % relative change of the solution
    ax(4,gapcounter) = nexttile;
    for m = 1:length(methods)
        if intervals
            h(m) = fillinterval(...
                2:maxit,...
                10.^means_4(gapcounter,2:maxit,m),...
                10.^lowers_4(gapcounter,2:maxit,m),...
                10.^uppers_4(gapcounter,2:maxit,m),...
                colors(m,:));
                hold on
        else
            h(m) = plot(1:maxit,10.^means_4(gapcounter,:,m),'color',colors(m,:));
            hold on
        end
    end
    legend(h,methodnames)
    xlabel('iteration')
    ylabel('relative solution change')
    set(gca,'YScale','log')
    title(sprintf('gap length %d ms',glengths(gapcounter)))
    xlim([0 50])
end
if intervals
    if bootstrap
        title(T,sprintf('bootstrap estimate of the mean performance from %d signals',length(signums)))
    else
        title(T,sprintf('normal estimate of the performance from %d signals',length(signums)))
    end
else
    title(T,sprintf('mean performance from %d signals',length(signums)))
end
for i = 1:4
    linkaxes(ax(i,:),'y')
end
set(gcf,'visible','on','windowstate','maximized')