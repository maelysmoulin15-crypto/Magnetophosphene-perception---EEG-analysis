function triggers_onset = detect_onsets_from_MF_auto( ...
    Lraw, fs, nb_expected, refractory_s, ignore_head_sec, env_smooth_sec, first_hint)
% detect_onsets_from_MF_auto
% Onsets via robust thresholds on a smoothed envelope and its derivative.
if nargin < 7, first_hint = []; end
if nargin < 6 || isempty(env_smooth_sec),  env_smooth_sec = 0.2; end
if nargin < 5 || isempty(ignore_head_sec), ignore_head_sec = 0.2; end
if nargin < 4 || isempty(refractory_s),    refractory_s = 8.0; end

x = double(Lraw(:));                         % column vector
N = numel(x);
w   = max(5, min(round(env_smooth_sec*fs), N));
env = smoothdata(abs(x), 'gaussian', w);     % smoothed |x|
denv = [0; diff(env)];                       % derivative

i0 = max(1, round(ignore_head_sec*fs)+1);    % ignore early clicks
seg_len = min(N - i0 + 1, round(10*fs));
if seg_len < 2000
    lo = max(1, i0-5000); hi = min(N, i0+5000); seg1 = env(lo:hi);
else
    seg1 = env(i0:i0+seg_len-1);
end
thr_env = median(seg1)      + 6*mad(seg1,1);
thr_dev = median(abs(denv)) + 6*mad(abs(denv),1);

cand = find( (env > thr_env) & (denv > thr_dev) & ((1:N)'>=i0) );

keep   = [];
last   = -inf;
minsep = round(refractory_s*fs);
for k = 1:numel(cand)
    if cand(k) - last >= minsep
        keep(end+1,1) = cand(k); %#ok<AGROW>
        last = cand(k);
    end
end

if ~isempty(first_hint)
    keep = unique([round(first_hint); keep], 'stable');
end
if ~isempty(nb_expected) && numel(keep) > nb_expected
    keep = keep(1:nb_expected);
end
triggers_onset = keep(:).';
end