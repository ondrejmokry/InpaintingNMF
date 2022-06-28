close all
clear
clc
addpath(genpath('utils'))

% suppress warnings
warning('off','all')

% file name
filename = 'results/inpainting_comparison_01.mat';
fprintf('Loading data...\n')
load(filename)

% prepare fields
PEMOQs    = cell(height(tables.(methods{1})),length(methods));
PEAQs     = cell(height(tables.(methods{1})),length(methods));
PEMOQgaps = NaN(height(tables.(methods{1})),1);
PEAQgaps  = NaN(height(tables.(methods{1})),1);

% folder with restored signals
sigfold = ['results/signals/',filename(9:end-4)];

%% processing
row = 0;
for signum = signums
    for glength = glengths

        row = row + 1;
        
        if ~isempty(PEAQs{row,length(methods)})
            continue
        end
        
        % command window output
        str = sprintf('Signal: %s',signals{signum});
        fprintf(repmat('=',length(str),1))
        fprintf('\n')
        fprintf(str)
        fprintf('\nGap length: %d ms\n',glength)
        fprintf(repmat('=',length(str),1))
        fprintf('\n')
              
        %% computing metrics
        for m = 1:length(methods)
            
            % loading the signals
            load(sprintf('%s/%s_%02d_%s.mat',sigfold,signals{signum}(1:3),glength,methods{m}),...
                'mask','signal','restored','fs')
            
            %% evaluate the gapped signal
            if m == 1
                % PEMO-Q
                fprintf('gapped signal, PEMO-Q\n')
                [~, ~, PEMOQgaps(row), ~] = audioqual(signal, signal.*mask, fs);
                fprintf(repmat('\b',1,22))
                
                % PEAQ
                fprintf('gapped signal, PEAQ\n')
                
                % save the reference signal as wav
                signal_48 = resample(signal, 48000, fs);
                audiowrite('signal_48.wav', signal_48, 48000);

                % save the gapped signal as wav
                gapped_48 = resample(signal.*mask, 48000, fs);
                audiowrite('gapped_48.wav', gapped_48, 48000);

                % evaluate
                PEAQgaps(row) = PQevalAudio_fn('signal_48.wav', 'gapped_48.wav', 0, length(gapped_48));
            end
            
            %% evaluate the reconstructed signal
            data_1 = NaN(size(restored,2),1);
            data_2 = NaN(size(restored,2),1);
            str = [];
            for i = 1:length(data_1)
                % PEMO-Q
                fprintf(repmat('\b',1,length(str)))
                str = sprintf('%s, iteration %d/%d\n',methods{m},i,length(data_1));
                fprintf(str)
                
                % evaluate
                [~, ~, data_1(i), ~] = audioqual(signal, restored(:,i), fs);
                fprintf(repmat('\b',1,22))
                
                % PEAQ                
                % save the restored signal as wav
                restored_48 = resample(restored(:,i), 48000, fs);
                audiowrite('restored_48.wav', restored_48, 48000);

                % evaluate
                data_2(i) = PQevalAudio_fn('signal_48.wav', 'restored_48.wav', 0, length(restored_48));
            end
            PEMOQs{row,m} = data_1;
            PEAQs{row,m}  = data_2;
            
        end

    end % gapnum

end % signum

%% save the data
for m = 1:length(methods)
    tables.(methods{m}).PEMOQgap = PEMOQgaps; 
    tables.(methods{m}).PEAQgap  = PEAQgaps;
    tables.(methods{m}).PEMOQ    = PEMOQs(:,m);
    tables.(methods{m}).PEAQ     = PEAQs(:,m);
end
save(filename,...
    'tables',...
    '-append')

%% clean up
delete gapped_48.wav restored_48.wav signal_48.wav