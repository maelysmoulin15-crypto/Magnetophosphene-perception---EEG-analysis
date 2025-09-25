function [triggers_onset, dbg] = detect_onsets_from_MF_auto( ...
    Lraw, fs, nb_expected, refractory_s, ignore_head_sec, env_smooth_sec)

% detect_onsets_from_MF_auto (100% automatique)
% - Upward threshold crossing on a smoothed |x| envelope
% - Must be followed by >= min_on_dur_s of sustained high envelope
% - Refractory period enforces at most 1 onset per exposure
% - If too few onsets, thresholds are relaxed iterativement
%
% INPUTS
%   Lraw            : raw detection signal (vector)
%   fs              : Hz
%   nb_expected     : expected number of real onsets (e.g. 50), [] if unknown
%   refractory_s    : min spacing (s) between two onsets (default 8)
%   ignore_head_sec : ignore first seconds (default 0.2)
%   env_smooth_sec  : envelope smoothing (s) (default 0.20)
%
% OUTPUTS
%   triggers_onset  : row vector of onset sample indices
%   dbg             : debug struct (envelope, thresholds, etc.)

    if nargin < 6 || isempty(env_smooth_sec),  env_smooth_sec  = 0.20; end
    if nargin < 5 || isempty(ignore_head_sec), ignore_head_sec = 0.20; end
    if nargin < 4 || isempty(refractory_s),    refractory_s    = 8.0;  end

    x = double(Lraw(:));
    N = numel(x);
    t = (0:N-1).'/fs;

    % --- envelope & derivative
    w   = max(5, min(round(env_smooth_sec*fs), N));
    env = smoothdata(abs(x), 'gaussian', w);
    denv = [0; diff(env)];

    % --- baseline (robuste): ne garder que les faibles valeurs d’enveloppe
    med_all = median(env);
    mad_all = mad(env,1);
    base_mask = env <= (med_all + 3*mad_all);
    base_env  = env(base_mask);
    if numel(base_env) < round(2*fs)   % fallback si trop court
        base_env = env;
    end

    % Seuils initiaux (serrés)
    K_env = 6; K_dev = 6;
    thr_env0 = median(base_env)      + K_env*mad(base_env,1);
    thr_dev0 = median(abs(denv))     + K_dev*mad(abs(denv),1);
    sustain_level = median(base_env) + 2*mad(base_env,1);  % niveau “haut”

    min_on_dur_s = 1.0;                              % doit rester haut ≥ 1 s
    sustain_need = round(0.8*min_on_dur_s*fs);       % tolérance 20 %
    refr_samp    = round(refractory_s*fs);
    i_start      = max(1, round(ignore_head_sec*fs)+1);

    % Boucle d’adaptation des seuils (si trop peu d’onsets)
    relax_steps = [1.0, 0.85, 0.7, 0.55, 0.45];      % on relâche progressivement
    best_idx = [];
    used_thr_env = thr_env0; used_thr_dev = thr_dev0;

    for r = 1:numel(relax_steps)
        thr_env = thr_env0 * relax_steps(r);
        thr_dev = thr_dev0 * relax_steps(r);

        % 1) franchissements vers le haut
        up = find(env(2:end) >= thr_env & env(1:end-1) < thr_env) + 1;
        up = up(up >= i_start);

        % 2) garder ceux suivis d’un “plateau” (env > sustain_level)
        keep = [];
        for k = 1:numel(up)
            i0 = up(k);
            i1 = min(N, i0 + round(2.5*fs));   % on vérifie sur ~2.5 s
            seg = env(i0:i1);
            if sum(seg > sustain_level) >= sustain_need
                keep(end+1,1) = i0; %#ok<AGROW>
            end
        end

        % 3) période réfractaire
        if ~isempty(keep)
            idx = keep(1);
            for k = 2:numel(keep)
                if keep(k) - idx(end) >= refr_samp
                    idx(end+1,1) = keep(k); %#ok<AGROW>
                end
            end
        else
            idx = [];
        end

        % 4) Si on a assez d’onsets, on s’arrête
        if isempty(nb_expected)
            best_idx = idx;
            used_thr_env = thr_env; used_thr_dev = thr_dev;
            break
        elseif numel(idx) >= nb_expected
            best_idx = idx(1:nb_expected);
            used_thr_env = thr_env; used_thr_dev = thr_dev;
            break
        else
            % sinon on retient le meilleur courant et on relaxe encore
            best_idx = idx;
            used_thr_env = thr_env; used_thr_dev = thr_dev;
        end
    end

    triggers_onset = best_idx(:).';   % row vector

    % Debug
    if nargout > 1
        dbg = struct('t',t,'env',env,'denv',denv, ...
                     'thr_env',used_thr_env,'thr_dev',used_thr_dev, ...
                     'sustain_level',sustain_level, ...
                     'refractory_s',refractory_s, ...
                     'min_on_dur_s',min_on_dur_s, ...
                     'onsets',triggers_onset);
    end
end
