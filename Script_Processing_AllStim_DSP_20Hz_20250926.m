%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                     EEG SIGNALS — MF EXPOSURE                       %%%
%%%                         — DSP ANALYSIS —                            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; clc; close all;

%% ----------------------- Paths / I/O ------------------------------------
addpath('/Users/maelys/Documents/MATLAB/');
addpath('/Users/maelys/Documents/MATLAB/fieldtrip');
ft_defaults;  % ajoute tous les sous-dossiers utiles

epoch_dir   = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/EEG_Magnetophosphes_Epoch';
results_dir = '/Users/maelys/Documents/Data/Magnétophosphènes/Results';
if ~exist(results_dir,'dir'), mkdir(results_dir); end
fig_dir = fullfile(results_dir, 'Figures_processing_PSD');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end

% (Optional) press table
press_xlsx  = '/Users/maelys/Documents/Data/Magnétophosphènes/Data/20Hz_global-CSV.xlsx';


%% ----------------------- Analysis parameters ----------------------------
bp_lo = 30;                      % Hz
bp_hi = 80;                      % Hz

notch_freqs   = [40 50 60];   % Hz
notch_bw_hz   = 3;               % bande totale (-2 dB) ~ ±0.5 Hz

t_keep        = [1 4];           % secondes conservées dans l’epoch 5 s

welch_win_sec = 1.0;             % s (mtmwelch + Hann)
welch_olap    = 0.5;             % (info) FieldTrip gère la segmentation interne

gamma_band    = [30 80];         % Hz
total_band    = [0 80];          % Hz (pour masques/figures)

% Donoghue aperiodic fit (inline)
fit_range     = [30 80];         % Hz, sans knee
exclude_bw_Hz = 1.0;             % exclure ±1 Hz autour des notches
z_thresh      = 2.5;             % sigma-clipping résidus positifs
max_iter      = 4;               % itérations robustes

% Figure settings
intensities_to_plot = 0:5:50;    % 11 colonnes
trial_of_interest   = 3;         % Trial 3
trial_tag = sprintf('_trial%d', trial_of_interest);
rows_per_fig        = 8;         % 8 électrodes par figure -> 8 figures / sujet
cols_per_fig        = numel(intensities_to_plot);

%% ----------------------- Progress settings ------------------------------
use_waitbar     = true;          % barre de progression
verbose         = true;          % messages console
progress_every  = 250;           % afficher toutes les N itérations
log_dir         = fullfile(results_dir, 'logs');
if ~exist(log_dir,'dir'), mkdir(log_dir); end
use_diary       = true;          % log texte
if use_diary
    diary(fullfile(log_dir, sprintf('dsp_log_%s.txt', datestr(now,'yyyymmdd_HHMMSS'))));
end

%% ----------------------- Load optional press table ----------------------
have_press = false;
if isfile(press_xlsx)
    Tpress = readtable(press_xlsx);
    Tpress.Properties.VariableNames = matlab.lang.makeValidName(Tpress.Properties.VariableNames);
    cname = Tpress.Properties.VariableNames;

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
        Tpress.(subj_col) = double(Tpress.(subj_col));
        Tpress.(inten_col)= double(Tpress.(inten_col));
        Tpress.(tr_col)   = double(Tpress.(tr_col));
    else
        warning('Press table présente mais colonnes non reconnues. Press = NaN.');
    end
end

%% ----------------------- Discover epoch files ---------------------------
L = dir(fullfile(epoch_dir, 'Sujet*_*Hz_*_*mT_trial*.mat'));
if isempty(L)
    error('Aucun fichier epoch .mat trouvé dans: %s', epoch_dir);
end
totalFiles = numel(L);
fprintf('Found %d epoch files\n', totalFiles);

% pré-scan sujets (pour logs)
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

rows = {};                 % accumulateur de table
psd_store = struct();      % pour figures

% --- init progression ---
t0 = tic;
if use_waitbar, hwb = waitbar(0, sprintf('Processing 0/%d files...', totalFiles)); end
last_subject = NaN;

%% ----------------------- Main loop over files ---------------------------
VarNames = {'Subject','Electrode','Intensity_mT','Trial', ...
    'Fs_original','Fs_effective', ...
    'GammaAbsPower_raw','GammaMeanFreq_raw', ...
    'GammaAbsPower_corr','GammaMeanFreq_corr', ...
    'AperiodicExponent','AperiodicOffset10', ...
    'EpochFile','Press'};

T = cell2table(cell(0, numel(VarNames)), 'VariableNames', VarNames);

for iFile = 1:10%totalFiles
    fpath = fullfile(L(iFile).folder, L(iFile).name);

    % Parse filename: Sujet%02d_%s_%dHz_%s_%dmT_trial%d.mat
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
    coil_exp = tok{2}; %#ok<NASGU>
    fq_stim  = str2double(tok{3}); %#ok<NASGU>
    electrode_label = tok{4};
    inten_mT = str2double(tok{5});
    trial_i  = str2double(tok{6});

    % --- log début de sujet ---
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

    % Charger fichier
    S = load(fpath);
    if ~isfield(S,'EEG_segment'), continue; end
    x_raw = double(S.EEG_segment(:));

    % Taux d’échantillonnage
    if isfield(S,'meta') && isfield(S.meta,'Fs')
        fs = double(S.meta.Fs);
    elseif isfield(S,'time_segment')
        tt = double(S.time_segment(:));
        fs = 1/median(diff(tt(~isnan(tt))));
    else
        continue;
    end
    if ~(isfinite(fs) && fs>0), continue; end

    % Retirer padding NaN éventuel
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


    %% ==== Filtrage Alex (Hann -> filtre 30–80) ====
% Hann sécurisé (colonnes) + option "periodic" (classique pour fenêtre d'analyse)
w        = hann(numel(x), 'periodic');
x_window = x(:) .* w;                 % x en colonne, même taille que w
x_epoch  = filtrealex(x_window, fs, 30, 0, 80, 0);   % -> signal déjà 30–80 Hz

%% ========= FieldTrip: ONLY notches (40/50/60) =========
data = [];
data.fsample = fs;
data.label   = {electrode_label};
data.trial   = {x_epoch(:)'};         % FT attend ligne
data.time    = {(0:numel(x_epoch)-1)/fs};

cfgp               = [];
cfgp.demean        = 'no';            % (ton band-pass Alex fixe déjà le DC)
% si tu préfères remettre à zéro la moyenne après Alex : 'yes'

cfgp.dftfilter     = 'yes';
cfgp.dftfreq       = [40 50 60];
cfgp.dftbandwidth  = notch_bw_hz;     % 1 Hz total
cfgp.dftreplace    = 'zero';

% !! pas de band-pass avec FT
cfgp.bpfilter      = 'no';

% padding inutile ici, on enlève
% cfgp.padding     = [];

data_prep = ft_preprocessing(cfgp, data);

% crop [1..4] s comme avant
cfgsel        = [];
cfgsel.latency= t_keep;
data_crop     = ft_selectdata(cfgsel, data_prep);


    %% ========= FOOOF aperiodic decomposition (robuste aux segments courts) =========
    % Segmenter façon Welch (2 s, 50% overlap) seulement si utile
    cfgseg = [];
    cfgseg.length  = 2;
    cfgseg.overlap = 0.5;

    % sécurité : ne segmenter que si on a assez d'échantillons pour 2 s
    seg_ok = false;
    try
        fs_loc = data_crop.fsample;
        nSamp  = numel(data_crop.time{1});
        seg_ok = isfinite(fs_loc) && fs_loc>0 && (nSamp/fs_loc) >= cfgseg.length;
    catch
        seg_ok = false;
    end

    if seg_ok
        data_seg = ft_redefinetrial(cfgseg, data_crop);
    else
        data_seg = data_crop; % pas assez long → pas de segmentation
    end

    % Durée effective (en s) du (premier) trial
    dur = numel(data_seg.time{1}) / data_seg.fsample;

    % Paramètres communs
    cfgS            = [];
    cfgS.method     = 'mtmfft';
    cfgS.foilim     = [30 80];
    cfgS.keeptrials = 'no';
    cfgS.pad        = 'nextpow2';

    % Choix du taper/smoothing selon la durée
    if dur >= 2.0
        % OK pour du multitaper lissé (DPSS)
        cfgS.taper     = 'dpss';
        cfgS.tapsmofrq = 1;     % lissage ~1 Hz (adapté à T>=2 s)
    elseif dur >= 1.0
        % Trop court pour DPSS lissé → FFT simple Hanning
        cfgS.taper     = 'hanning';
        if isfield(cfgS,'tapsmofrq'), cfgS = rmfield(cfgS,'tapsmofrq'); end
    else
        % Vraiment trop court → ignorer ce fichier (évite les erreurs FOOOF/PSD)
        warning('Trial trop court (%.3fs) → ignoré: %s', dur, L(iFile).name);
        continue;
    end

    % Apériodique via FOOOF
    cfgS.output = 'fooof_aperiodic';
    fractal     = ft_freqanalysis(cfgS, data_seg);

    % Spectre total (mêmes paramètres)
    cfgS.output = 'pow';
    original    = ft_freqanalysis(cfgS, data_seg);

    % Composante oscillatoire = total - apériodique
    cfgm = [];
    cfgm.parameter = 'powspctrm';
    cfgm.operation = 'x2-x1';
    oscillatory    = ft_math(cfgm, fractal, original);

    % --- métriques 30–80 (sur la grille d’original) ---
   % --- métriques 30–80 (sur la grille d’original) ---
f    = original.freq(:);
Pap  = squeeze(fractal.powspctrm);   Pap  = Pap(:);
Ptot = squeeze(original.powspctrm);  Ptot = Ptot(:);

% masque gamma
mask_gamma = (f>=30 & f<=80);

% masque autour des notches (même largeur que le notch DFT)
bw_tot   = notch_bw_hz;          % p.ex. 3 Hz TOTAL
half_bw  = bw_tot/2;
notch_mask = false(size(f));
for fc = notch_freqs(:).'
    notch_mask = notch_mask | (f > (fc-half_bw) & f < (fc+half_bw));
end

% composante oscillatoire = total - apériodique
Posc = Ptot - Pap;

% sécurité : pas d’énergie négative après soustraction
Posc(Posc < 0) = 0;

% bins valides pour métriques/figures
valid = mask_gamma & ~notch_mask;

df = mean(diff(f));
gamma_abs_raw    = sum(Ptot(valid))              * df;
gamma_fmean_raw  = sum(f(valid).*Ptot(valid))    * df / max(gamma_abs_raw,eps);
gamma_abs_corr   = sum(Posc(valid))              * df;
gamma_fmean_corr = sum(f(valid).*Posc(valid))    * df / max(gamma_abs_corr,eps);

% ---- Fit du 1/f sur la partie apériodique (log-log), en excluant notches
fit_valid = valid & Pap>0 & isfinite(Pap);
if any(fit_valid)
    p = polyfit(log10(f(fit_valid)), log10(Pap(fit_valid)), 1); % LP = p(1)*LF + p(2)
    aperiodic_exponent = -p(1);        % convention Donoghue (slope négative -> exponent positive)
    aperiodic_offset10 = 10.^p(2);     % offset à 10^intercept
else
    aperiodic_exponent = NaN;
    aperiodic_offset10 = NaN;
end

% pour les figures
f_band        = f(valid);
Pxx_band_raw  = Ptot(valid);
Pxx_band_corr = Posc(valid);




    % Optional press
    press_val = NaN;
    if have_press
        rows_press = Tpress.(subj_col)==subj & ...
            Tpress.(inten_col)==inten_mT & ...
            Tpress.(tr_col)==trial_i;
        if any(rows_press)
            press_val = Tpress.(press_col)(find(rows_press,1));
        end
    end

    % -------- Accumulate row --------
    rowT = cell2table({ ...
        subj, electrode_label, inten_mT, trial_i, fs, fs, ...
        gamma_abs_raw, gamma_fmean_raw, ...
        gamma_abs_corr, gamma_fmean_corr, ...
        aperiodic_exponent, aperiodic_offset10, ...
        L(iFile).name, press_val}, ...
        'VariableNames', VarNames);
    T = [T; rowT]; %#ok<AGROW>

    % -------- Stockage pour figures --------
    if contains(L(iFile).name, trial_tag)
        f_band        = f(mask_gamma);
        Pxx_band_raw  = Ptot(mask_gamma);
        Pxx_band_corr = Posc(mask_gamma);

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

    % -------- Progression --------
    if use_waitbar && (mod(iFile,progress_every)==0 || iFile==totalFiles)
        try
            waitbar(iFile/totalFiles, hwb, ...
                sprintf('Processing %d/%d...', iFile, totalFiles));
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
if isempty(fieldnames(psd_store))
    warning('Aucune donnée stockée dans psd_store → pas de figures générées.');
else
    subjects = fieldnames(psd_store);
    for s = 1:numel(subjects)
        Skey = subjects{s};
        elec_names = sort(fieldnames(psd_store.(Skey)));

        % Détecter toutes les intensités disponibles chez ce sujet
        all_I_keys = {};
        for e = 1:numel(elec_names)
            Ekey = elec_names{e};
            all_I_keys = [all_I_keys; fieldnames(psd_store.(Skey).(Ekey))]; %#ok<AGROW>
        end
        all_I_keys = unique(all_I_keys);
        I_all = sort(str2double(erase(all_I_keys, 'I')));

        % Choix des colonnes : intersection avec ta liste cible sinon toutes
        I_target = intersect(intensities_to_plot, I_all);
        if isempty(I_target), I_target = I_all; end
        cols_per_fig = max(1, numel(I_target));

        nElec = numel(elec_names);
        nFigs = ceil(nElec / rows_per_fig);
        fprintf('Composing figures for %s: %d electrodes -> %d figures (intensités: %s)\n', ...
            Skey, nElec, nFigs, mat2str(I_target));

        for figi = 1:nFigs
            rStart = (figi-1)*rows_per_fig + 1;
            rEnd   = min(figi*rows_per_fig, nElec);

            fig = figure('Name', sprintf('%s — PSD 30–80 (trial %d)', Skey, trial_of_interest), ...
                'Units','pixels','Position',[50 50 2000 1200], 'Visible','off');
            tl = tiledlayout(fig, rows_per_fig, cols_per_fig, 'Padding','compact','TileSpacing','compact');
            sgtitle(tl, sprintf('%s — PSD 30–80 Hz (raw & 1/f-corrected) — trial %d', ...
                Skey, trial_of_interest));

            for rr = 1:rows_per_fig
                eIdx = rStart + rr - 1;
                for cc = 1:cols_per_fig
                    inten = I_target(cc);
                    Ikey  = sprintf('I%03d', inten); % clé robuste avec 3 chiffres
                    nexttile(tl, (rr-1)*cols_per_fig + cc);

                    if eIdx <= nElec
                        Ekey = elec_names{eIdx};
                        hasData = isfield(psd_store.(Skey).(Ekey), Ikey);
                        if hasData
                            D = psd_store.(Skey).(Ekey).(Ikey);
                            f_plot = D.f;
                            y1 = 10*log10(D.raw + eps);
                            y2 = 10*log10(D.corr + eps);
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

            % Sauvegarde
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

