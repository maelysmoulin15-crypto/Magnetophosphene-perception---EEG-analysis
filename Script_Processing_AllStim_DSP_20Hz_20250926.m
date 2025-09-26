%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     EEG SIGNALS — MF EXPOSURE                       %%%
%%%                         — DSP ANALYSIS —                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all;   % Clean workspace, clear command window, close figures

%% ----------------------- Paths / I/O ------------------------------------
addpath('/Users/maelys/Documents/MATLAB/');                 % Path that contains 'filtrealex'
addpath('/Users/maelys/Documents/MATLAB/fieldtrip');        % FieldTrip root folder
ft_defaults;                                                % Initialize FieldTrip paths and defaults

epoch_dir   = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/EEG_Magnetophosphes_Epoch'; % Folder with per-epoch .mat files
results_dir = '/Users/maelys/Documents/Data/Magnétophosphènes/Results';                         % Output folder for tables and .mat results
fig_dir = fullfile(results_dir, 'Figures_processing_PSD');                                      % Subfolder for PSD figures

press_xlsx  = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/20Hz_global-CSV.xlsx';      % (Optional) behavioral/press table


%% ----------------------- Analysis parameters ----------------------------
bp_lo = 30;                      % Lower bound of the (conceptual) gamma band of interest [Hz]
bp_hi = 80;                      % Upper bound of the gamma band of interest [Hz]

notch_freqs   = [40 50 60];      % Power-line and harmonics to be suppressed [Hz]
notch_bw_hz   = 3;               % Display/semantic bandwidth reference (not used directly later) [Hz]

t_keep        = [1 4];           % Time window (in seconds) extracted from each 5-s epoch for analysis

welch_win_sec = 1.0;             % For PSD estimation, an informative reference; actual segmentation handled by FT
welch_olap    = 0.5;             % Fractional overlap (informative; FT manages internally)

gamma_band    = [30 80];         % Band used for metrics (area, mean frequency) [Hz]
total_band    = [0 80];          % Global band for masking/figures [Hz]

% Donoghue (aperiodic 1/f) fit — we request FT's "fooof_aperiodic" output later
fit_range     = [30 80];         % Fitting range in Hz (no knee model here)
exclude_bw_Hz = 1.0;             % To exclude ±1 Hz around notch frequencies when fitting/plotting (handled explicitly later)
z_thresh      = 2.5;             % Robust residual clipping threshold (informative; not used below)
max_iter      = 4;               % Max iterations for robust fitting (informative; not used below)

% Figure settings
intensities_to_plot = 0:5:50;    % Columns (stimulus intensities) included in summary figures
trial_of_interest   = 3;         % Trial index to visualize
trial_tag = sprintf('_trial%d', trial_of_interest);  % String pattern for picking trials for figures
rows_per_fig        = 8;         % Number of electrodes per figure row block
cols_per_fig        = numel(intensities_to_plot);    % Number of columns equals number of intensities plotted

%% ----------------------- Progress settings ------------------------------
use_waitbar     = true;          % GUI progress bar
verbose         = true;          % Console progress messages
progress_every  = 250;           % Print a message every N files
log_dir         = fullfile(results_dir, 'logs');  % Text logs folder
if ~exist(log_dir,'dir'), mkdir(log_dir); end
use_diary       = true;          % Whether to record console output to a log file
if use_diary
    diary(fullfile(log_dir, sprintf('dsp_log_%s.txt', datestr(now,'yyyymmdd_HHMMSS'))));
end

%% ----------------------- Load optional press table ----------------------
% Tries to read an optional behavioral table. If present and recognizable,
% its values can be joined into the final results table.
have_press = false;
if isfile(press_xlsx)
    Tpress = readtable(press_xlsx);
    Tpress.Properties.VariableNames = matlab.lang.makeValidName(Tpress.Properties.VariableNames);
    cname = Tpress.Properties.VariableNames;

    % Resolve column names that may vary (French/English)
    if ismember('Sujet',cname),       subj_col = 'Sujet';
    elseif ismember('Subject',cname), subj_col = 'Subject';
    else, subj_col = ''; end

    if ismember('intensite',cname),       inten_col = 'intensite';
    elseif ismember('Intensity',cname),   inten_col = 'Intensity';
    else, inten_col = ''; end

    tr_col = '';
    for k = 1:numel(cname)
        if any(strcmpi(cname{k},{'essai','Trial','Essai'})), tr_col = cname{k}; break; end
    end

    press_col = '';
    for k = 1:numel(cname)
        if any(strcmpi(cname{k},{'press','Press','PRESS','Pression'})), press_col = cname{k}; break; end
    end

    have_press = ~isempty(subj_col) && ~isempty(inten_col) && ~isempty(tr_col) && ~isempty(press_col);
    if have_press
        % Ensure numeric types for matching
        Tpress.(subj_col) = double(Tpress.(subj_col));
        Tpress.(inten_col)= double(Tpress.(inten_col));
        Tpress.(tr_col)   = double(Tpress.(tr_col));
    else
        warning('Press table present but required columns not recognized. Will use NaN for press.');
    end
end

%% ----------------------- Discover epoch files ---------------------------
% Expect filenames like: Sujet%02d_%s_%dHz_%s_%dmT_trial%d.mat
L = dir(fullfile(epoch_dir, 'Sujet*_*Hz_*_*mT_trial*.mat'));
if isempty(L)
    error('No epoch .mat files found in: %s', epoch_dir);
end
totalFiles = numel(L);
fprintf('Found %d epoch files\n', totalFiles);

% Pre-scan to count files per subject (for nicer logs)
subjects_all = nan(totalFiles,1);
for iFile = 1:totalFiles
    tok = regexp(L(iFile).name, '^Sujet(\d+)_', 'tokens','once');
    if ~isempty(tok), subjects_all(iFile) = str2double(tok{1}); end
end
subj_unique = unique(subjects_all(~isnan(subjects_all)));
subj_counts = containers.Map('KeyType','double','ValueType','double');
for s = 1:numel(subj_unique)
    ss = subj_unique(s);
    subj_counts(ss) = sum(subjects_all==ss);
end

rows = {};            % (Unused accumulator placeholder; kept for compatibility)
psd_store = struct(); % Storage for spectra used in summary figures

% Progress handle and bookkeeping
t0 = tic;
if use_waitbar, hwb = waitbar(0, sprintf('Processing 0/%d files...', totalFiles)); end
last_subject = NaN;

%% ----------------------- Main loop over files ---------------------------
% Final results table columns
VarNames = {'Subject','Electrode','Intensity_mT','Trial', ...
    'Fs_original','Fs_effective', ...
    'GammaAbsPower_raw','GammaMeanFreq_raw', ...
    'GammaAbsPower_corr','GammaMeanFreq_corr', ...
    'AperiodicExponent','AperiodicOffset10', ...
    'EpochFile','Press'};

T = cell2table(cell(0, numel(VarNames)), 'VariableNames', VarNames);

for iFile = 1:10%totalFiles    % NOTE: limited to first 10 for testing; set to totalFiles to process everything
    fpath = fullfile(L(iFile).folder, L(iFile).name);

    % Parse filename pieces (subject, coil, stim frequency, electrode label, intensity, trial index)
    tok = regexp(L(iFile).name, ...
        '^Sujet(\d+)_([A-Za-z]+)_(\d+)Hz_(.+)_(\d+)mT_trial(\d+)\.mat$', ...
        'tokens','once');
    if isempty(tok)
        if verbose && mod(iFile,progress_every)==0
            fprintf('[%6d/%6d] Skip (bad name): %s\n', ...
                iFile, totalFiles, L(iFile).name);
        end
        continue;
    end
    subj     = str2double(tok{1});
    coil_exp = tok{2}; %#ok<NASGU>   % Coil/experiment tag not used directly here
    fq_stim  = str2double(tok{3}); %#ok<NASGU>   % Stimulation frequency (not used directly)
    electrode_label = tok{4};
    inten_mT = str2double(tok{5});
    trial_i  = str2double(tok{6});

    % Pretty logging when a new subject starts
    if ~isequal(subj, last_subject)
        if verbose
            total_this_subj = 0;
            if isKey(subj_counts, subj)
                total_this_subj = subj_counts(subj);
            end
            fprintf('\n=== Subject %02d: %d files to process ===\n', ...
                subj, total_this_subj);
        end
        last_subject = subj;
    end

    % ----------------------- Load epoch -----------------------
    S = load(fpath);
    if ~isfield(S,'EEG_segment'), continue; end
    x_raw = double(S.EEG_segment(:));   % Single-channel epoch

    % Sampling rate inference (prefer meta.Fs; else derive from time vector)
    if isfield(S,'meta') && isfield(S.meta,'Fs')
        fs = double(S.meta.Fs);
    elseif isfield(S,'time_segment')
        tt = double(S.time_segment(:));
        fs = 1/median(diff(tt(~isnan(tt))));
    else
        continue;
    end
    if ~(isfinite(fs) && fs>0), continue; end

    % Remove any padding NaNs (epochs may be stored with tail NaNs)
    if isfield(S,'meta') && isfield(S.meta,'short_epoch_pad') && ...
            ~isempty(S.meta.short_epoch_pad)
        goodN = numel(x_raw) - double(S.meta.short_epoch_pad);
        goodN = max(0, min(numel(x_raw), goodN));
        last_valid = max(goodN, find(~isnan(x_raw),1,'last'));
    else
        tmp_last = find(~isnan(x_raw),1,'last');
        if isempty(tmp_last), continue; end
        last_valid = tmp_last;
    end
    x = x_raw(1:last_valid);

    %% ==== Alex band-pass stage (user-defined) ============================
    % Apply a Hann window prior to the custom filter, as required by 'filtrealex'.
    % IMPORTANT: This is the *only* band-pass stage (30–80 Hz). FieldTrip is later
    % used *only* for notches and spectral analysis.
    w        = hann(numel(x), 'periodic');     % Analysis-friendly Hann
    x_window = x(:) .* w;                      % Column vector multiplication
    x_epoch  = filtrealex(x_window, fs, 30, 0, 80, 0);   % Custom band-pass to 30–80 Hz

    %% ========= FieldTrip: NOTCH ONLY (40/50/60 Hz) ======================
    % Build a minimal FT raw structure to run filters on the already band-passed data.
    data = [];
    data.fsample = fs;
    data.label   = {electrode_label};
    data.trial   = {x_epoch(:)'};              % FT expects 1xN row vectors for trials
    data.time    = {(0:numel(x_epoch)-1)/fs};

    % We use windowed-sinc FIR band-stop filters (stable, linear phase).
    % The order scales with sampling rate and half-bandwidth to give strong attenuation.
    bw  = 1;  % half-width in Hz → e.g., 50±1 Hz stopband   (the next token keeps the original text)
    bw  = 0.5;  % demi-largeur (→ stop 50±1 Hz)             % <-- kept exactly as in your script

    cfgp = [];
    cfgp.demean        = 'no';                    % DC already handled by Alex stage; keep as-is
    cfgp.bsfilter      = 'yes';                   % Enable band-stop
    cfgp.bsfilttype    = 'firws';                 % Windowed-sinc FIR
    cfgp.bsfreq        = [40-bw 40+bw; 50-bw 50+bw; 60-bw 60+bw];   % Three notches
    cfgp.bsfiltwintype = 'hamming';               % Hamming window for FIR design
    cfgp.bsfiltord     = max(ceil(3*fs/bw), 500); % Heuristic: higher order → deeper notch
    cfgp.bpfilter      = 'no';                    % No additional FT band-pass
    data_notch = ft_preprocessing(cfgp, data);    % Apply the notches

    % ----------------------- Crop to analysis window ----------------------
    % Keep only the central portion [1..4] s to avoid transient edges.
    cfgsel        = [];
    cfgsel.latency= t_keep;
    data_crop     = ft_selectdata(cfgsel, data_notch);

    %% ========= Spectral analysis & aperiodic (1/f) via FT ===============
    % Optionally segment (Welch-like) if long enough (2 s, 50% overlap).
    cfgseg = [];
    cfgseg.length  = 2;
    cfgseg.overlap = 0.5;

    % Only segment if we have enough samples for a 2-s window
    seg_ok = false;
    try
        fs_loc = data_crop.fsample;
        nSamp  = numel(data_crop.time{1});
        seg_ok = isfinite(fs_loc) && fs_loc>0 && (nSamp/fs_loc) >= cfgseg.length;
    catch
        seg_ok = false;
    end

    if seg_ok
        data_seg = ft_redefinetrial(cfgseg, data_crop); % Multiple overlapped trials
    else
        data_seg = data_crop; % Use single trial without segmentation
    end

    % Effective duration (first trial); used to choose taper below
    dur = numel(data_seg.time{1}) / data_seg.fsample;

    % Common spectral config (multitaper if long enough, else Hann)
    cfgS            = [];
    cfgS.method     = 'mtmfft';
    cfgS.foilim     = [30 80];     % Restrict to gamma range
    cfgS.keeptrials = 'no';
    cfgS.pad        = 'nextpow2';  % Zero-padding to next power of 2 for finer frequency grid

    if dur >= 2.0
        cfgS.taper     = 'dpss';   % Multitaper with ~1 Hz smoothing
        cfgS.tapsmofrq = 1;
    elseif dur >= 1.0
        cfgS.taper     = 'hanning'; % Simpler taper when segments are short
        if isfield(cfgS,'tapsmofrq'), cfgS = rmfield(cfgS,'tapsmofrq'); end
    else
        % If we have <1 s, spectra and FOOOF are unreliable → skip this file
        warning('Trial too short (%.3fs) → skipped: %s', dur, L(iFile).name);
        continue;
    end

    % Ask FieldTrip to return the aperiodic component estimated by FOOOF
    cfgS.output = 'fooof_aperiodic';
    fractal     = ft_freqanalysis(cfgS, data_seg);

    % Compute the total power spectrum with the same settings
    cfgS.output = 'pow';
    original    = ft_freqanalysis(cfgS, data_seg);

    % Oscillatory component = total - aperiodic
    cfgm = [];
    cfgm.parameter = 'powspctrm';
    cfgm.operation = 'x2-x1';          % original - fractal
    oscillatory    = ft_math(cfgm, fractal, original);

    % === Metrics & figure vectors in 30–80 Hz (excluding notched bins) ===
    f    = original.freq(:);
    Ptot = squeeze(original.powspctrm); Ptot = Ptot(:);
    Pap  = squeeze(fractal.powspctrm);  Pap  = Pap(:);
    Posc = squeeze(oscillatory.powspctrm); Posc = Posc(:);

    % Build a frequency mask for the gamma band and exclude exact notch regions
    mask_gamma = (f>=30 & f<=80);
    bw = 2;  % Total width for exclusion: [39–41], [49–51], [59–61]
    notch_mask = (abs(f-40)<bw/2) | (abs(f-50)<bw/2) | (abs(f-60)<bw/2);
    valid = mask_gamma & ~notch_mask;

    % Safety: clamp tiny negative oscillatory values to 0 (numerical noise)
    Posc(Posc<0) = 0;

    % Area (absolute power) and mean frequency within the valid gamma bins
    df               = mean(diff(f(valid)));
    gamma_abs_raw    = sum(Ptot(valid))*df;
    gamma_fmean_raw  = sum(f(valid).*Ptot(valid))*df / max(gamma_abs_raw,eps);
    gamma_abs_corr   = sum(Posc(valid))*df;
    gamma_fmean_corr = sum(f(valid).*Posc(valid))*df / max(gamma_abs_corr,eps);

    % Recover classic 1/f parameters from the aperiodic curve (log10–log10 linear fit)
    aperiodic_exponent = NaN;
    aperiodic_offset10 = NaN;
    vv = valid & Pap>0 & isfinite(Pap);
    if any(vv)
        p = polyfit(log10(f(vv)), log10(Pap(vv)), 1);  % log10(P) = p1*log10(f) + p2
        aperiodic_exponent = -p(1);                   % by convention exponent = -slope
        aperiodic_offset10 = 10.^p(2);               % offset at 10^0 = 1 Hz in linear units
    end

    % Store vectors (already excluding notched bins) for plotting later
    f_band        = f(valid);
    Pxx_band_raw  = Ptot(valid);
    Pxx_band_corr = Posc(valid);

    % -------- Optional press (behavioral) value join --------
    press_val = NaN;
    if have_press
        rows_press = Tpress.(subj_col)==subj & ...
            Tpress.(inten_col)==inten_mT & ...
            Tpress.(tr_col)==trial_i;
        if any(rows_press)
            press_val = Tpress.(press_col)(find(rows_press,1));
        end
    end

    % -------- Append a row to the results table --------
    rowT = cell2table({ ...
        subj, electrode_label, inten_mT, trial_i, fs, fs, ...
        gamma_abs_raw, gamma_fmean_raw, ...
        gamma_abs_corr, gamma_fmean_corr, ...
        aperiodic_exponent, aperiodic_offset10, ...
        L(iFile).name, press_val}, ...
        'VariableNames', VarNames);
    T = [T; rowT]; %#ok<AGROW>

    % -------- Save spectra for figures (selected trial only) --------
    if contains(L(iFile).name, trial_tag)
        Skey = sprintf('S%02d', subj);
        Ekey = matlab.lang.makeValidName(electrode_label);
        Ikey = sprintf('I%03d', inten_mT);
        if ~isfield(psd_store, Skey), psd_store.(Skey) = struct(); end
        if ~isfield(psd_store.(Skey), Ekey), psd_store.(Skey).(Ekey) = struct(); end
        psd_store.(Skey).(Ekey).(Ikey) = struct( ...
            'f',   f_band, ...
            'raw', Pxx_band_raw, ...
            'corr',Pxx_band_corr );
    end

    % -------- Progress feedback --------
    if use_waitbar && (mod(iFile,progress_every)==0 || iFile==totalFiles)
        try
            waitbar(iFile/totalFiles, hwb, sprintf('Processing %d/%d...', iFile, totalFiles));
        end
    end
    if verbose && (mod(iFile,progress_every)==0 || iFile==1 || iFile==totalFiles)
        elapsed = toc(t0);
        rate = iFile / max(elapsed, eps);
        remain = (totalFiles - iFile) / max(rate, eps);
        fprintf('[%6d/%6d] %5.1f%% | elapsed %6.1fs | remain %6.1fs | subj %02d\n', ...
            iFile, totalFiles, 100*iFile/totalFiles, elapsed, remain, subj);
        drawnow limitrate;
    end
end

%% ----------------------- Build and save table ---------------------------
% Sort and write results; also save a .mat copy for convenient re-use.
if height(T)==0
    warning('No results collected. Check your epoch directory and file naming.');
    if use_diary, diary off; end
    return;
end

T = sortrows(T, {'Subject','Electrode','Intensity_mT','Trial'});

csv_out = fullfile(results_dir, 'DSP_Results_Gamma_30_80_withDonoghue.csv');
mat_out = fullfile(results_dir, 'DSP_Results_Gamma_30_80_withDonoghue.mat');
writetable(T, csv_out);
save(mat_out, 'T');

fprintf('\nDone metrics. Wrote:\n  %s\n  %s\n', csv_out, mat_out);

%% ----------------------- Compose & save the figures ---------------------
% Compose tiled figures: rows are electrodes, columns are intensities.
if isempty(fieldnames(psd_store))
    warning('No spectra stored in psd_store → no figures produced.');
else
    subjects = fieldnames(psd_store);
    for s = 1:numel(subjects)
        Skey = subjects{s};
        elec_names = sort(fieldnames(psd_store.(Skey)));

        % Discover all intensities available for this subject (any electrode)
        all_I_keys = {};
        for e = 1:numel(elec_names)
            Ekey = elec_names{e};
            all_I_keys = [all_I_keys; fieldnames(psd_store.(Skey).(Ekey))]; %#ok<AGROW>
        end
        all_I_keys = unique(all_I_keys);
        I_all = sort(str2double(erase(all_I_keys, 'I')));

        % Use requested list if present; otherwise plot all
        I_target = intersect(intensities_to_plot, I_all);
        if isempty(I_target), I_target = I_all; end
        cols_per_fig = max(1, numel(I_target));

        nElec = numel(elec_names);
        nFigs = ceil(nElec / rows_per_fig);
        fprintf('Composing figures for %s: %d electrodes -> %d figures (intensities: %s)\n', ...
            Skey, nElec, nFigs, mat2str(I_target));

        for figi = 1:nFigs
            rStart = (figi-1)*rows_per_fig + 1;
            rEnd   = min(figi*rows_per_fig, nElec);

            % Hidden figure for batch export
            fig = figure('Name', sprintf('%s — PSD 30–80 (trial %d)', Skey, trial_of_interest), ...
                'Units','pixels','Position',[50 50 2000 1200], 'Visible','off');
            tl = tiledlayout(fig, rows_per_fig, cols_per_fig, 'Padding','compact','TileSpacing','compact');
            sgtitle(tl, sprintf('%s — PSD 30–80 Hz (raw & 1/f-corrected) — trial %d', ...
                Skey, trial_of_interest));

            for rr = 1:rows_per_fig
                eIdx = rStart + rr - 1;
                for cc = 1:cols_per_fig
                    inten = I_target(cc);
                    Ikey  = sprintf('I%03d', inten);
                    nexttile(tl, (rr-1)*cols_per_fig + cc);

                    if eIdx <= nElec
                        Ekey = elec_names{eIdx};
                        hasData = isfield(psd_store.(Skey).(Ekey), Ikey);
                        if hasData
                            D = psd_store.(Skey).(Ekey).(Ikey);
                            f_plot = D.f;
                            % Plot in linear units (µV^2/Hz). Switch to dB with 10*log10 if desired.
                            plot(f_plot, D.raw); hold on; plot(f_plot, D.corr);
                            ylabel('\muV^2/Hz');
                            xlim([30 80]); grid on;
                        else
                            text(0.5,0.5,'NA','HorizontalAlignment','center','FontSize',10);
                            axis off;
                        end
                    else
                        axis off;
                    end

                    if rr==1, title(sprintf('%dmT', inten), 'FontSize',9); end
                    if cc==1 && eIdx <= nElec
                        ylabel(strrep(elec_names{eIdx},'_','\_'), 'Interpreter','tex', 'FontSize',9);
                    end
                    if rr==1 && cc==1
                        legend({'Raw','Corrected (1/f)'}, 'Location','southoutside', ...
                            'Orientation','horizontal', 'Box','off');
                    end
                    set(gca,'FontSize',8);
                end
            end

            % Export both PNG and MATLAB FIG for later edits
            png_out = fullfile(fig_dir, sprintf('%s_PSD30-80_trial%d_part%d.png', ...
                Skey, trial_of_interest, figi));
            fig_out = fullfile(fig_dir, sprintf('%s_PSD30-80_trial%d_part%d.fig', ...
                Skey, trial_of_interest, figi));
            exportgraphics(fig, png_out, 'Resolution', 200);
            savefig(fig, fig_out);
            close(fig);
            fprintf('  Saved %s\n', png_out);
            drawnow limitrate;
        end
    end
    fprintf('\nFigures saved in: %s\n', fig_dir);
end
