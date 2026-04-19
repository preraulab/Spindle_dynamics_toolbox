function quick_start
%QUICK_START  Launch a GUI to specify a spindle-dynamics model and generate an overview figure
%
%   Usage:
%       quick_start
%
%   Inputs:
%       none
%
%   Outputs:
%       none (side effects only - opens a uifigure and produces an overview figure)
%
%   Notes:
%       GUI steps:
%         1) load data (bundled example or user data)
%         2) choose model factors (SOphase, stage, SOpower, history)
%         3) choose history order if history is selected
%         4) choose interaction terms
%       Clicking "Run the Model" triggers the pipeline and renders the
%       overview figure.
%
%       To load user data, the .mat file must contain EEG, Fs, stage_val,
%       and stage_time (variable names are matched case-insensitively):
%         EEG        : 'EEG', 'eeg_data', 'raw_EEG', 'data'
%         Fs         : 'Fs', 'sampling_rate'
%         stage_val  : 'stage_val', 'stage_vals', 'stages', 'stage', 'stageval', 'stagevals'
%         stage_time : 'stage_time', 'stage_t', 'time_stages', 'stage_times'
%
%       Accompanies Chen et al., PNAS 2025.
%
%   See also: Example_Script, specify_mdl, preprocessToDesignMatrix
%
%   ∿∿∿  Prerau Laboratory MATLAB Codebase · sleepEEG.org  ∿∿∿

%% GUI   
    % Add path
    addpath(genpath('./helper_function'))

    % Initialize global cache for res_table so you dont run it everytime
    res_table = [];
    
    % Create a separate floating window for options
    optionsFig = uifigure('Name', 'Quick Start GUI', 'Position', [1050, 100, 300, 500]);

    % 1. Load Data:
    uilabel(optionsFig, 'Position', [20, 450, 200, 20], 'Text', '1. Load Data:');
    btnExampleData = uibutton(optionsFig, 'Text', 'Load Example Data', ...
        'Position', [20, 420, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) loadExampleData());
    btnUserData = uibutton(optionsFig, 'Text', 'Load User Data', ...
        'Position', [150, 420, 120, 30], ...
        'ButtonPushedFcn', @(btn, event) loadUserData());

    % Display data status
    dataStatusLabel = uilabel(optionsFig, 'Position', [20, 390, 260, 20], ...
        'Text', 'No data loaded', 'HorizontalAlignment', 'left');

    % 2. Checkboxes for selecting factors
    uilabel(optionsFig, 'Position', [20, 360, 200, 20], 'Text', '2. Select Factors:');
    cbSOphase = uicheckbox(optionsFig, 'Text', 'SOphase', 'Position', [20, 330, 100, 20]);
    cbStage = uicheckbox(optionsFig, 'Text', 'Stage', 'Position', [120, 330, 100, 20]);
    cbSOpower = uicheckbox(optionsFig, 'Text', 'SOpower', 'Position', [20, 300, 100, 20]);
    cbHistory = uicheckbox(optionsFig, 'Text', 'History', 'Position', [120, 300, 100, 20], ...
        'ValueChangedFcn', @(src, event) toggleHistoryChoice());

    % 3. Radio buttons for history choice
    uilabel(optionsFig, 'Position', [20, 270, 200, 20], 'Text', '3. Select History Order:');
    bgHist = uibuttongroup(optionsFig, 'Position', [20, 220, 200, 40], 'Enable', 'off'); % Initially disabled
    rbShort = uiradiobutton(bgHist, 'Text', 'Short-term (15 sec)', 'Position', [10, 20, 150, 20]);
    rbLong = uiradiobutton(bgHist, 'Text', 'Long-term (90 sec)', 'Position', [10, 0, 150, 20]);
    rbLong.Value = true; % Default to "long"

    % 4. Dropdown for interaction terms
    uilabel(optionsFig, 'Position', [20, 190, 200, 20], 'Text', '4. Select Interactions:');
    ddInteract = uidropdown(optionsFig, 'Position', [20, 160, 200, 22], ...
        'Items', {'None', 'stage:SOphase', 'SOpower:SOphase', 'stage:SOphase, stage:history'}, ...
        'Value', 'None'); % Default to "None"

    % 5. "Run Model" button
    btnRun = uibutton(optionsFig, 'Text', 'Run the Model', ...
                      'Position', [20, 40, 150, 30], ...
                      'ButtonPushedFcn', @(btn, event) runModel());

    % Initialize variables to store data
    dataLoaded = false;
    EEG = [];
    Fs = [];
    stage_val = [];
    stage_time = [];

    % Function to toggle History Choice availability
    function toggleHistoryChoice()
        if cbHistory.Value
            bgHist.Enable = 'on';  % Enable history choice options
        else
            bgHist.Enable = 'off'; % Disable history choice options
        end
    end

    % Function to load example data
    function loadExampleData()
        data_path = 'example_data/';
        load([data_path, 'example_data1.mat'], 'EEG', 'Fs', 'stage_val', 'stage_time'); % Replace with your actual variables
        res_table = []; % Reset cached res_table since new data is loaded
        dataLoaded = true;
        dataStatusLabel.Text = 'Example data loaded successfully!';
    end

    
    % Function to load user data
    function loadUserData()
        [file, path] = uigetfile('*.mat', 'Select Data File');
        if isequal(file, 0)
            dataStatusLabel.Text = 'User canceled data selection.';
        else
            try
                loadedData = load(fullfile(path, file)); % Load user-selected file
                
                % Check for potential variable names for each required field
                EEG_names = {'EEG', 'eeg_data', 'raw_EEG','data'}; % Possible names for EEG
                Fs_names = {'Fs', 'sampling_rate'}; % Possible names for Fs
                stage_val_names = {'stage_val', 'stage_vals', 'stages','stage','stageval','stagevals'}; % Possible names for stage values
                stage_time_names = {'stage_time', 'stage_t', 'time_stages','stage_t','stage_times'}; % Possible names for stage times
    
                % Initialize variables
                EEG = [];
                Fs = [];
                stage_val = [];
                stage_time = [];
    
                % Check each field in loadedData
                fields = fieldnames(loadedData);
    
                % Match and assign the appropriate variable
                for i = 1:length(fields)
                    field = fields{i};
                    if any(strcmpi(field, EEG_names))
                        EEG = double(loadedData.(field));
                    elseif any(strcmpi(field, Fs_names))
                        Fs = double(loadedData.(field));
                    elseif any(strcmpi(field, stage_val_names))
                        stage_val = double(loadedData.(field));
                    elseif any(strcmpi(field, stage_time_names))
                        stage_time = double(loadedData.(field));
                    end
                end
    
                % Verify that all required variables are assigned
                if isempty(EEG) || isempty(Fs) || isempty(stage_val) || isempty(stage_time)
                    uialert(optionsFig, 'Missing required variables: EEG, Fs, stage_val, or stage_time.', 'Load Error');
                else
                    % If all variables are assigned successfully
                    dataLoaded = true;
                    res_table = []; % Reset cached res_table since new data is loaded
                    dataStatusLabel.Text = ['User data loaded: ', file];
                end
    
            catch ME
                uialert(optionsFig, ['Error loading file: ', ME.message], 'Load Error');
            end
        end
    end

    % Function to run the model and generate a separate figure
    function runModel()
        % Ensure data is loaded
        if ~dataLoaded
            uialert(optionsFig, 'No data loaded. Please load example or user data.', 'Data Missing');
            return;
        end

        % Get user inputs
        factors = [cbSOphase.Value, cbStage.Value, cbSOpower.Value, cbHistory.Value];
        selectedInteraction = ddInteract.Value;
        switch selectedInteraction
            case 'None'
                InteractSelect = {};
            case 'stage:SOphase'
                InteractSelect = {'stage:SOphase'};
            case 'SOpower:SOphase'
                InteractSelect = {'SOpower:SOphase'};
            case 'stage:SOphase, stage:history'
                InteractSelect = {'stage:SOphase', 'stage:history'};
        end
        hist_choice = 'long'; % Default to 'long'
        if rbShort.Value
            hist_choice = 'short';
        end

        % Validate inputs
        if ~any(factors)
            uialert(optionsFig, 'Please select at least one factor.', 'Input Error');
            return;
        end

        % Check if `res_table` is cached
        if isempty(res_table)
            disp('Extracting spindle info ...');
            res_table = eegToEventSignal(EEG, Fs, stage_val, stage_time); % Compute only once
        end

        % Call preprocessing and fitting function
        BinarySelect = double(factors); % Convert logical to numeric
        [X, BinData, ModelSpec] = preprocess_GUI(BinarySelect, InteractSelect, hist_choice, EEG, Fs, stage_val, stage_time,res_table);

        % Generate the overview figure in a new window
        plot_overview_GUI(EEG, res_table, X, BinData, ModelSpec);
    end
end



% ****************************HELPER FUNCTIONS*************************
% Preprocess Function
function [X, BinData, ModelSpec] = preprocess_GUI(BinarySelect, InteractSelect, hist_choice, EEG, Fs, stage_val, stage_time,res_table)
    % Ensure required inputs are provided
    if isempty(EEG) || isempty(Fs) || isempty(stage_val) || isempty(stage_time)
        error('EEG, Fs, stage_val, and stage_time must be provided.');
    end

    disp('Data preprocessing started... ');

    % Preprocess EEG data
    binsize = 0.1; % Bin size in seconds
    hard_cutoffs = [12, 16]; % Frequency range for fast spindles

    % Compute spindle events and their properties only when it is empty
    if isempty(res_table)
        disp('May take minutes ... Extracting spindles ...');
        res_table = eegToEventSignal(EEG, Fs, stage_val, stage_time);
    end
    
    % Convert event data to binned data
    [BinData] = rawToBinData(res_table, Fs, binsize, hard_cutoffs, hist_choice);
    disp('Data preprocessing complete.');

    % Define control points based on history choice
    switch lower(hist_choice)
        case 'short'
            hist_ord = 15 / binsize;
            control_pt = [0:15:90 120 hist_ord];
        case 'long'
            hist_ord = 90 / binsize;
            control_pt = [0:15:90 120 150:100:750 hist_ord];
    end

    % Specify model based on user inputs
    ModelSpec = specify_mdl(BinarySelect, InteractSelect, ...
                             'hist_choice', hist_choice, ...
                             'control_pt', control_pt);

    % Build design matrix for the specified model
    X = build_design_mt(ModelSpec, BinData);
    disp('Design Matrix complete.');
    
end


% Plot Overview Function
function plot_overview_GUI(EEG,res_table,X, BinData, ModelSpec)
    % Prepare data for figure
    binsize = ModelSpec.binsize;
    xtime = BinData.real_time + binsize/2;        % bin data time stamps 
    sta = BinData.sta;                            % stage binary matrix [N1,N2,N3,REM,Wake]  
    stagetrain = 3*sta(:,1)+ 2*sta(:,2)+ sta(:,3)+4*sta(:,4)+5*sta(:,5); % stagetrain in bin
    y = BinData.y;                                % spindle binary train
    sop = BinData.sop;                            % SO power in bin  
    phase = BinData.phase;                        % SO phase in bin  
    Fs = BinData.Fs;
    hard_cutoffs = ModelSpec.hard_cutoffs;

    % Compute multitaper spectrogram
    dsfreqs = 0.1;
    nfft = 2^(nextpow2(Fs/dsfreqs));
    [mt, stimes, sfreqs] = multitaper_spectrogram_mex(EEG, Fs,[8,19],[2,3],[1,0.05],nfft,'constant',[],false,false);
    
    % Extract fast spindle raw time and freq
    fast_inds = (res_table.peak_freqs{:} >= hard_cutoffs(1)) & (res_table.peak_freqs{:}< hard_cutoffs(2));
    spindletime_all = res_table.peak_ctimes{:};
    freq_all  = res_table.peak_freqs{:};
    fast_pt = spindletime_all(fast_inds);
    fast_pf = freq_all(fast_inds);
    
    % Extract raw SO power and phase
    phase_raw = cell2mat(res_table.SOphase);      % SO phase raw
    sop_raw = cell2mat(res_table.SOpow);          % SO power raw

    % Fit GLM Model
    disp('Model fitting ...')
    [b, ~, stats] = glmfit(X, y, 'poisson');
   
    %***************** Generate the overview figure*************
    disp('Figure generating ...')
    figure('Name', 'Spindle Dynamics Overview', 'Position', [100, 50, 1400, 1000]);

    % --------------------------- Left Panel ---------------
    % 1) Spectrogram
    ax1 = subplot(7, 3, [1 2 4 5]); % Merge subplots
    hold on;
    imagesc(stimes, sfreqs, nanpow2db(mt));
    axis xy;
    ylabel('Frequency (Hz)');
    climscale;
    c = colorbar_noresize;
    colormap(gca, rainbow4);
    ylabel(c, 'Power (dB)');
    axis tight;
    plot(fast_pt, fast_pf, 'o', 'markersize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    title('Sleep EEG Spectrogram');
    scaleline(ax1, 10, '10 s');
    xlabel(" "+newline+" ")
    
    % 2) Spindle train
    ax2 = subplot(7, 3, [7 8]); % Merge subplots
    stem(fast_pt, ones(size(fast_pt)), 'Color', 'k', 'Marker', 'none', 'LineWidth', 1.5);
    set(ax2,'XTick', [],'YColor', 'none','box','off');
    ylim(ax2, [0 1.3]);
    title('Sleep Spindle Train');
    xlabel(" "+newline+"   ")

    % 3) Stage train
    ax3 = subplot(7, 3, [10 11]); % Merge subplots
    hold on;
    hypnoplot(xtime, stagetrain);
    plot(xtime(y == 1), stagetrain(y == 1), 'o', 'markersize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    title('Sleep Stage');
    xlabel(" "+newline+"   ")
    set(ax3, 'XTick', [])

    % 4) SO power
    ax5 = subplot(7, 3, [13 14]); % Merge subplots
    hold on;
    t = (0:length(EEG)-1)/Fs;
    plot(t, sop_raw, 'LineWidth', 1.5, 'Color', 'k');
    plot(xtime(y == 1), sop(y == 1), 'o', 'markersize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    ylabel('Power');
    title('Continuous Sleep Depth (SO Power)');
    set(ax5, 'XTick', [])

    % 5) SO phase
    ax6 = subplot(7, 3, [16 17]); % Merge subplots
    hold on;
    plot(t, phase_raw, 'k', 'LineWidth', 1.5);
    phase = wrapToPi(phase);
    plot(xtime(y == 1), phase(y == 1), 'o', 'markersize', 6, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
    set(ax6, 'XTick', [],'YTick', [-pi 0 pi], 'YTickLabels', {'-\pi', '0', '\pi'}, 'ylim', [-pi pi]);
    yline(0, 'k--');
    ylabel('Phase');
    title('Cortical Up/Down State (SO Phase)');
   

    % 6) CIF
    ax7 = subplot(7, 3, [19 20]); % Merge subplots
    [CIF_sta, ~, ~] = glmval(b, X, 'log', stats); % Compute instantaneous spindle rate
    hold on;
    plot(xtime, CIF_sta * 60 / binsize, 'color', 'k', 'LineWidth', 2);
    plot(xtime(y == 1), CIF_sta(y == 1) * 60 / binsize, 'o', 'markersize', 6, ...
         'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'k', 'LineWidth', 2);
    title(sprintf('Density = f (%s)', ModelSpec.FactorSelect));
    xlabel('Time (s)', 'FontSize', 12);
    ylabel({'Spindle Density'; '(events/min)'}, 'FontSize', 12);
    
    % Adjust position for CIF plot
    pos = get(ax7, 'Position');
    pos(4) = pos(4) + 0.05;
    pos(2) = pos(2) - 0.05;
    set(ax7, 'Position', pos);
    %scaleline(ax7, 5, '5 secs');
    
    % Link x-axes and scrollzoompan
    linkaxes([ax1, ax2, ax3, ax5, ax6, ax7], 'x');
    scrollzoompan(ax7, 'x');

    %--------------------- Right panel -----------------------
    % 1) Plot Polarhistogram
    ax4 = subplot(7, 3, [3 6 9]); % Define subplot location
    phasehistogram(phase(y==1),1,'NumBins',50); 
    rlim([0 .4])
    rticks(0:0.1:0.4)
    title('Spindle SO Phase Histogram',FontSize=12)
    
    % 2) Plot history Curve
    ax8 = subplot(7, 3, [12 15 18]); % Define subplot location
    if ModelSpec.BinarySelect(4)==1
        plot_hist_curve(stats,ModelSpec,BinData);
    end
        
    % Store axes handles in a cell array
    axesHandles = {ax1, ax2, ax3, ax5, ax6, ax7};
    
    % Loop through each axis and adjust properties
    for i = 1:length(axesHandles)
        currentAx = axesHandles{i};
        title(currentAx, currentAx.Title.String, 'FontSize', 10, 'FontWeight', 'bold'); 
        set(currentAx, 'FontSize', 10, 'box', 'off');
    end

end
