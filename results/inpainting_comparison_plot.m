clear
clc
close all
addpath('../utils')

% settings
filename  = 'inpainting_comparison_01.mat';
bootstrap = false;
intervals = false;

% load data
fprintf('Loading %s...\n',filename)
load(filename)
signals = signals(signums);

% reorder according to the submitted paper
methods = methods([1, 2, 3, 8, 9, 7, 5, 4, 6, 10]);
methodnames = methodnames([1, 2, 3, 8, 9, 7, 5, 4, 6, 10]);

% measures
measures = {'SNR','PEMOQ','PEAQ'};
ylabels  = {'SNR (dB)','PEMO-Q ODG','PEAQ ODG'};
variants = {'end','max'};

% prepare figure for SNR and ODG
f_data = figure('visible','off','windowstate','maximized');
T_data = tiledlayout(length(variants),length(measures),'tileindexing','columnmajor');
colors = colororder();
if length(methods) > 7
    M = length(methods)-7;
    colors = [colors; [zeros(M,2), (1:M)'/M]];
end
colororder(colors)

% prepare figure for time
f_time = figure('visible','off','windowstate','maximized');
T_time = tiledlayout(length(variants),length(measures),'tileindexing','columnmajor');
colororder(colors)

for measure = 1:length(measures)
    
    %% reorganize the data for plotting
    str = sprintf('measure: %s',measures{measure});
    fprintf(repmat('=',1,length(str)))
    fprintf('\n%s\n',str)
    fprintf(repmat('=',1,length(str)))
    fprintf('\n')
    fprintf('Reorganizing the data...\n')
    DATend = NaN(length(signums),length(glengths),length(methods));
    DATmax = NaN(length(signums),length(glengths),length(methods));
    TIMend = NaN(length(signums),length(glengths),length(methods));
    TIMmax = NaN(length(signums),length(glengths),length(methods));
    for i = 1:length(signums)
        for j = 1:length(glengths)
            for k = 1:length(methods)

                % find the row
                rows = strcmp(tables.(methods{k}).signal,signals{i});
                rows = rows .* (tables.(methods{k}).gap == glengths(j));
                row  = find(rows);
                if isempty(row)
                    continue
                end

                % read the value and save
                DATend(i,j,k) = tables.(methods{k}).(measures{measure}){row}(end);
                [DATmax(i,j,k), I] = max(tables.(methods{k}).(measures{measure}){row});
                TIMend(i,j,k) = tables.(methods{k}).time(row);
                TIMmax(i,j,k) = tables.(methods{k}).time(row) * I/length(tables.(methods{k}).(measures{measure}){row});
            end
        end
    end

    %% plotting means
    for fig = 1:2
        ax = gobjects(length(variants),length(measures));
        for v = 1:length(variants)
            if fig == 1
                ax(v,measure) = nexttile(T_data);
            else
                ax(v,measure) = nexttile(T_time);
            end
            fprintf('data: %s, variant: %s\n',measures{measure},variants{v})
            if strcmpi(variants{v},'end')
                if fig == 2
                    data = TIMend;
                else
                    data = DATend;
                end
            else
                if fig == 2
                    data = TIMmax;
                else
                    data = DATmax;
                end
            end

            if intervals
                if bootstrap
                    % interval estimate using bootstrapping
                    fprintf('   Computing the bootstraps...\n') %#ok<*UNRCH>
                    [means, lowers, uppers] = bootstrap_est(data);
                else
                    fprintf('   Computing the stds...\n')
                    mult   = tinv(0.975,length(signums)-1);
                    means  = squeeze(mean(data,1,'omitnan'));
                    stds   = squeeze(std(data,0,1,'omitnan'));
                    lowers = means - mult*stds/sqrt(length(signums));
                    uppers = means + mult*stds/sqrt(length(signums));
                end
            else
                means = squeeze(mean(data,1,'omitnan'));
            end

            % plot
            if intervals
                fprintf('   Plotting...\n')
            end
            h = gobjects(length(methods),1);
            for m = 1:length(methods)
                if intervals
                    h(m) = fillinterval(glengths,means(:,m),lowers(:,m),uppers(:,m),colors(m,:));
                else
                    h(m) = plot(glengths,means(:,m));
                    hold on
                end
            end
            legend(h,methods)
            xlabel('gap length (ms)')
            if fig == 2
                set(gca,'YScale','log')
                ylabel('elapsed time (s)')
            else
                ylabel(ylabels{measure})
            end
            if fig == 2
                title(sprintf('elapsed time for %s, variant: %s',measures{measure},variants{v}))
            else
                title(sprintf('%s, variant: %s',measures{measure},variants{v}))
            end
            ylim('padded')
        end
        linkaxes(ax(:,measure),'y')
    end
    set(f_data,'visible','on')
    set(f_time,'visible','on')

    %% plotting individual signals
    fprintf('Plotting individual results...\n')
    for v = 1:length(variants)
        f = figure('visible','off','windowstate','maximized');
        colororder(colors)
        t = tiledlayout(2,5,'tileindexing','columnmajor');
        ax = gobjects(length(signums),1);

        for s = 1:length(signums)
            ax(s) = nexttile;
            if strcmpi(variants{v},'end')
                data = DATend;
            else
                data = DATmax;
            end

            % plot
            plot(glengths,squeeze(data(s,:,:)))
            legend(methodnames)
            xlabel('gap length (ms)')
            ylabel(ylabels{measure})
            title(signals{s},'interpreter','none')

        end

        % link y axis and show
        linkaxes(ax,'y')
        title(t,sprintf('%s, variant: %s',measures{measure},variants{v}))
        set(f,'visible','on')
    end
end