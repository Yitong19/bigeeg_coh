%% This script is for calculating single parameters for behavioral data and rejecting erroneous trials
INDIR = 'C:\Users\yitong_z\Documents\BIGEEG\data\beh_preprocess';
mat_files = dir(fullfile(INDIR, '*.mat'));
beh_names = {mat_files.name};

for isubj = 1:length(beh_names)

    path_subj = fullfile(INDIR,beh_names{isubj});
    A = load(path_subj).beh_clean;

    for b = 1:length(A)

        n_block = A{1,b};

        for ti = 1:length(n_block)

            current_trial = n_block{1,ti};

            if size(current_trial,1) > 0

                % frequency by fast Fourrier transform
                y = double(current_trial.y);
                y = y - mean(y);  % Remove the constant offset
                N = numel(y);
                
                Y = fft(y);

                % One-sided amplitude spectrum
                amplitude = abs(Y(1:floor(N/2)+1)) / N;

                if rem(N, 2) == 0
                    amplitude(2:end-1) = 2 * amplitude(2:end-1);
                else
                    amplitude(2:end) = 2 * amplitude(2:end);
                end
                
                % Corresponding frequencies
                f = (0:floor(N/2))' * Fs / N;
                
                figure;
                plot(f, amplitude);
                xlabel('Frequency (Hz)');
                ylabel('Amplitude');
                xlim([0 10]);
                grid on;



            else

            end


        end


    end






end
































