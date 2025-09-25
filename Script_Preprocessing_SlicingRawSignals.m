%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     EEG SIGNALS — MF EXPOSURE                       %%%
%%%                       — SLICING RAW SIGNALS —                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script:
%   1) Scans all subjects on disk (pattern "SujetXX_<EXP>_20Hz.cnt").
%   2) Loads the order of intensities (55 trials) from a .txt (1 column/subject).
%   3) Auto-selects a DETECTION channel (priority to M1/TP9/A1/LM… else MF/Laser).
%   4) Detects onsets of real exposures on that RAW detection channel.
%   5) Builds the 55 epochs of 5 s while STRICTLY respecting the .txt order:
%        - real (≠0 mT):    start = onset(k) of the k-th real
%        - sham (0 mT):     5-s window placed ±10 s around a real onset (no overlap)
%   6) Saves the 55 RAW epochs for all EEG electrodes (64 if MF=65).
%   7) Exports 1 figure “RAW detection channel” + 4 figures of 16 EEG subplots.
%
% Notes:
%   • No filtering or analysis is performed here (raw segments are saved).
%   • Epoch length policy at file edges:
%       - epochs are ONSET-ALIGNED;
%       - if the tail is missing, it is padded (NaN by default) to 5 s.
%
% PREREQUISITES
%   - MATLAB + EEGLAB (for pop_loadcnt / Neuroscan .cnt)
%   - Order file "RandomSubject_<EXP>.txt" (size 55 x N_subjects)
%   - Utility functions at the bottom of this file or as separate .m files:
%       functions/pick_detection_channel.m
%       functions/detect_onsets_from_MF_auto.m
%       functions/build_epochs_from_order.m
%
% AUTHOR  : Maëlys MOULIN
% VERSION : 1.1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all;

%% ========================= PATHS / SETTINGS ========================== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% File naming parameters (as in the .cnt filenames)
fq_stim  = 20;           % e.g., "20Hz"
coil_exp = "GLO";        % e.g., "GLO" (use the same token as in filenames)

% Add your toolboxes
addpath('/Users/maelys/Documents/MATLAB/');
addpath('/Users/maelys/Documents/MATLAB/eeglab2025.0.0');
eeglab;  % load I/O (no need to open the UI)

% I/O locations
raw_data_path   = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/EEG_Magnetophosphes_RawData';
stim_order_file = sprintf('/Users/maelys/Documents/Data/Magnétophosphènes/Data/RandomSubject_%s.txt', coil_exp);

seg_base_dir = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/EEG_Magnetophosphes_Epoch';
fig_base_dir = '/Users/maelys/Documents/Data/Magnétophosphènes/Results';
fig_dir      = fullfile(fig_base_dir, 'Figures_preprocessing');

if ~exist(seg_base_dir,'dir'), mkdir(seg_base_dir); end
if ~exist(fig_base_dir,'dir'), mkdir(fig_base_dir); end
if ~exist(fig_dir,'dir'),     mkdir(fig_dir);      end

%% ============== DETECTION / EPOCH / FIGURE PARAMETERS =============== %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Onset detection (on the RAW detection channel)
refractory_s     = 8.0;    % minimum spacing between onsets (s) ~ inter-stim
ignore_head_sec  = 0.2;    % ignore first 200 ms (start clicks)
env_smooth_sec   = 0.20;   % envelope |x| smoothing (s) for detection

% Epoch construction
epoch_sec        = 5;      % epoch duration (s)
sham_margin_sec  = 10;     % sham windows placed ±10 s around a real onset

% Epoch completion policy at file edges (keeps onset alignment)
pad_short_epochs   = true; % pad tail if not enough samples left
pad_value          = NaN;  % use NaN (or 0) for padding
require_full_epoch = false;% set true to drop short epochs instead of padding

% Figures
make_laser_fig   = true;   % detection channel + onsets + epoch markers
make_grid16      = true;   % 4 figures of 16 EEG subplots (up to 64 channels)
grid_ds_maxpts   = 5000;   % ~points per channel for display (downsampling)

% (OPTION) If the very first onset must be imposed for some subjects:
% Example: manual_first_onset.S35 = 14.3172 % seconds (Subject 35)
manual_first_onset = struct();  % leave empty for 100% auto
% manual_first_onset.S35 = 14.327;
% manual_first_onset.S41 = 2.4728;
% manual_first_onset.S44 = 5.2338;


%% ================== LOAD ORDER FILE + LIST SUBJECTS ================= %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfile(stim_order_file)
    error('Order file not found: %s', stim_order_file);
end
stim_orders = readmatrix(stim_order_file);  % expected 55 x N_subjects

% List available .cnt
pat = sprintf('Sujet*_%s_%dHz.cnt', coil_exp, fq_stim);
L = dir(fullfile(raw_data_path, pat));
if isempty(L)
    error('No files "%s" in %s', pat, raw_data_path);
end

% Extract subject numbers
subs = [];
rx = sprintf('^Sujet(\\d+)_%s_\\d+Hz\\.cnt$', coil_exp);  % anchor on current EXP
for k = 1:numel(L)
    tok = regexp(L(k).name, rx, 'tokens','once');
    if ~isempty(tok), subs(end+1) = str2double(tok{1}); end %#ok<SAGROW>
end
subs = unique(sort(subs));
fprintf('=== %d subject(s) found: %s ===\n', numel(subs), num2str(subs));

% subs currently contains all subjects found on disk
only_these = [35 41 44];               % <- choose who to rerun
subs = intersect(subs, only_these);    % keep only those that exist on disk
fprintf('=== Re-running subjects: %s ===\n', num2str(subs));

%% ============================ MAIN LOOP ============================= %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s = 1:numel(subs)
    sub = subs(s);

    % -------- Load EEG --------
    filename = sprintf('Sujet%02d_%s_%dHz.cnt', sub, coil_exp, fq_stim);
    file = fullfile(raw_data_path, filename);
    if ~isfile(file)
        fprintf('Subject %02d: missing file -> SKIP\n', sub);
        continue;
    end

    fprintf('\n>>> Subject %02d : loading %s\n', sub, filename);
    eeg = pop_loadcnt(file, 'dataformat','auto');
    fs  = eeg.srate;
    N   = size(eeg.data,2);
    time = (0:N-1)/fs;

    % -------- Get this subject’s order (55 trials) --------
    if sub > size(stim_orders,2)
        warning('Subject %02d: no order column in %s -> SKIP', sub, stim_order_file);
        continue;
    end
    intensities_order = stim_orders(:, sub);
    if numel(intensities_order) ~= 55
        error('Subject %02d: order column must contain 55 values.', sub);
    end
    nb_onset_expected = sum(intensities_order ~= 0); % typically 50

    % -------- Choose detection channel --------
    det_idx   = pick_detection_channel(eeg);
    det_label = eeg.chanlocs(det_idx).labels;
    fprintf('Detection channel: %s (idx %d)\n', det_label, det_idx);

    % -------- EEG channels to slice (exclude MF if present) --------
    nCh      = size(eeg.data,1);
    labelsUp = upper(string({eeg.chanlocs.labels}));
    mf_keys  = ["LASER","MF","MAG","COIL","STIM","EXPOSURE"];
    mf_found = find(ismember(labelsUp, mf_keys), 1);

    if nCh >= 65 && ~isempty(mf_found) && mf_found == 65
        eeg_cols = 1:64;                            % 64 EEG + 65=MF
    else
        eeg_cols = setdiff(1:nCh, mf_found);        % all except MF if identified
    end

    % -------- Detect onsets on RAW detection signal (100% auto) --------
    Det = double(eeg.data(det_idx,:));
    Det = detrend(Det);
    Det = Det - mean(Det);

    % --- Ignorer un préfixe uniquement pour certains sujets (détection seule)
    ignore_head_override = struct('S35', 14.0);  % secondes à ignorer pour S35
    key = sprintf('S%02d', sub);
    if isfield(ignore_head_override, key)
        ignore_head_local = ignore_head_override.(key);
        fprintf('   Detection: ignoring first %.3f s for subject %02d\n', ignore_head_local, sub);
    else
        ignore_head_local = ignore_head_sec;  % valeur globale (0.2 s)
    end

    [triggers_onset, dbg] = detect_onsets_from_MF_auto( ...
        Det, fs, nb_onset_expected, refractory_s, ignore_head_local, env_smooth_sec);

    % -------- Build 55 epoch indices (respect .txt order) --------
    [i_Beg_stim, i_End_stim] = build_epochs_from_order( ...
        triggers_onset, intensities_order, fs, N, epoch_sec, sham_margin_sec);
    % NOTE: i_End_stim is clamped to N (can be short at the end by design)

    % -------- Figure: RAW detection channel + markers --------
    if make_laser_fig
        ds = max(1, floor(N/20000));              % ~20k points for display
        tt = time(1:ds:end);

        fig = figure('Name', sprintf('Subject %02d — %s RAW (detection)', sub, det_label), ...
            'Units','pixels','Position',[60 80 1600 500],'Visible','on');
        plot(time, Det, 'k'); hold on; grid on; axis tight;
        yl = ylim;

        % onsets = red triangles
        yv = yl(1) + 0.95*(yl(2)-yl(1));
        plot(triggers_onset/fs, yv*ones(size(triggers_onset)), 'rv', ...
            'MarkerFaceColor','r', 'MarkerSize',5);

        % vertical lines = epoch starts (colored by intensity)
        t_sham   = time(i_Beg_stim(intensities_order == 0));
        t_50     = time(i_Beg_stim(intensities_order == 50));
        t_others = time(i_Beg_stim(~ismember(intensities_order,[0 50])));

        line([t_sham; t_sham],   repmat(yl(:),1,numel(t_sham)),   'LineStyle','--','Color',[1 0 1],   'LineWidth',0.8);
        line([t_50;   t_50],     repmat(yl(:),1,numel(t_50)),     'LineStyle','--','Color',[1 0 0],   'LineWidth',0.8);
        line([t_others; t_others],repmat(yl(:),1,numel(t_others)),'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.5);


        legend({'Raw detection channel','Onsets','0 mT','50 mT','Others'}, ...
            'Location','northoutside','Orientation','horizontal');
        xlabel('Time (s)'); ylabel('Amplitude (raw)');
        title(sprintf('Subject %02d — detection channel: %s', sub, det_label));

        pngL = sprintf('Fig_DetectionRaw_Subject%02d_%s_%dHz.png', sub, coil_exp, fq_stim);
        figL = sprintf('Fig_DetectionRaw_Subject%02d_%s_%dHz.fig', sub, coil_exp, fq_stim);
        exportgraphics(fig, fullfile(fig_dir, pngL), 'Resolution', 300);
        savefig(fig, fullfile(fig_dir, figL));
        close(fig);
    end

    % -------- Save 55 RAW epochs for ALL EEG channels --------
    Lepoch_samp = round(epoch_sec * fs);

    for channel = eeg_cols
        % Channel label (fallback to "ChXX") and a safe filename token
        lab = eeg.chanlocs(channel).labels;
        if isempty(lab), lab = sprintf('Ch%02d', channel); end
        lab_safe = regexprep(lab, '[^\w-]+', '');

        x_raw = double(eeg.data(channel,:));

        for i = 1:55
            i0 = i_Beg_stim(i);
            i1 = i_End_stim(i);

            if i0 > N || i1 < 1, continue; end  % out of bounds (very unlikely)

            % extract (may be shorter than Lepoch_samp at the end of file)
            i0c = max(1, i0);
            i1c = min(N, i1);
            seg = x_raw(i0c:i1c);
            tsg = time(i0c:i1c);

            % pad tail if needed (keeps onset alignment)
            Lseg = numel(seg);
            if Lseg < Lepoch_samp
                if require_full_epoch
                    % skip short epochs
                    fprintf('   Skip short epoch: ch %d, trial %d (%d < %d samples)\n', ...
                        channel, i, Lseg, Lepoch_samp);
                    continue;
                elseif pad_short_epochs
                    tmp = repmat(pad_value, Lepoch_samp, 1);
                    tmp(1:Lseg) = seg;
                    seg = tmp;

                    ttmp = NaN(Lepoch_samp, 1);
                    ttmp(1:Lseg) = tsg;
                    tsg = ttmp;
                end
            end

            EEG_segment  = seg;
            time_segment = tsg;
            inten        = intensities_order(i);

            % Trial index 1..5 within THIS intensity (per channel)
            trial_num = sum(intensities_order(1:i) == inten);

            % Meta
            meta = struct( ...
                'subject', sub, ...
                'fq_stim', fq_stim, ...
                'electrode_label', lab, ...
                'electrode_index', channel, ...
                'Fs', fs, ...
                'epoch_sec', epoch_sec, ...
                'intensity_mT', inten, ...
                'trial_num', trial_num, ...     % 1..5 within this intensity
                'trial_abs', i, ...             % absolute order 1..55
                'det_channel', det_label, ...
                'short_epoch_pad', max(0, Lepoch_samp - Lseg), ...
                'pad_value', pad_short_epochs * pad_value ...
                );

            seg_name = sprintf('Sujet%02d_%s_%dHz_%s_%dmT_trial%d.mat', ...
                sub, coil_exp, fq_stim, lab_safe, inten, trial_num);
            filepath = fullfile(seg_base_dir, seg_name);
            save(filepath, 'EEG_segment', 'time_segment', 'inten', 'trial_num', 'meta');
        end
    end

    % -------- 4 figures of 16 EEG subplots (with markers) --------
    if make_grid16
        chans = eeg_cols(1:min(64, numel(eeg_cols)));
        nper  = 16;
        nfig  = ceil(numel(chans)/nper);
        ds    = max(1, floor(N/max(1,grid_ds_maxpts)));
        tt    = time(1:ds:end);

        t_sham   = time(i_Beg_stim(intensities_order == 0));
        t_50     = time(i_Beg_stim(intensities_order == 50));
        t_others = time(i_Beg_stim(~ismember(intensities_order,[0 50])));

        for fi = 1:nfig
            idx = (fi-1)*nper + 1 : min(fi*nper, numel(chans));
            fig = figure('Name', sprintf('Subject %02d — channels %d–%d', sub, idx(1), idx(end)), ...
                'Units','pixels','Position',[60 60 1600 900],'Visible','on');
            tl = tiledlayout(fig, 4, 4, 'Padding','compact','TileSpacing','compact');
            sgtitle(tl, sprintf('Subject %02d — Figure %d/%d', sub, fi, nfig));

            for k = 1:numel(idx)
                ch  = chans(idx(k));
                lab = eeg.chanlocs(ch).labels; if isempty(lab), lab = sprintf('Ch%02d', ch); end
                x   = double(eeg.data(ch, 1:ds:end));

                nexttile;
                plot(tt, x, 'k'); hold on; grid on; axis tight;
                yl = ylim;
                line([t_sham; t_sham],   repmat(yl(:),1,numel(t_sham)),   'LineStyle','--','Color',[1 0 1],   'LineWidth',0.8);
                line([t_50;   t_50],     repmat(yl(:),1,numel(t_50)),     'LineStyle','--','Color',[1 0 0],   'LineWidth',0.8);
                line([t_others; t_others],repmat(yl(:),1,numel(t_others)),'LineStyle','--','Color',[0.7 0.7 0.7],'LineWidth',0.5);

                title(lab, 'Interpreter','none', 'FontSize',8);
                set(gca,'XTick',[],'YTick',[]);
            end

            pngG = sprintf('Fig_AllElectrodsDetection_Subject%02d_%s_%dHz_part%d.png', sub, coil_exp, fq_stim, fi);
            figG = sprintf('Fig_AllElectrodsDetection_Subject%02d_%s_%dHz_part%d.fig', sub, coil_exp, fq_stim, fi);
            exportgraphics(fig, fullfile(fig_dir, pngG), 'Resolution', 300);
            savefig(fig, fullfile(fig_dir, figG));
            close(fig);
        end
    end

    fprintf('>>> OK: Subject %02d — %d EEG channels processed (detection: %s)\n', ...
        sub, numel(eeg_cols), det_label);
end

fprintf('\n=== Done ===\n');


