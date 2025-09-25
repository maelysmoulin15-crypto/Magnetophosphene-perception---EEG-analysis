function [i_Beg_stim, i_End_stim] = build_epochs_from_order( ...
    triggers_onset, intensities_order, fs, N, epoch_sec, sham_margin_sec)
% build_epochs_from_order
% Returns start/end indices for the 55 epochs (5 s), respecting the .txt order.
Lepoch = round(epoch_sec * fs);
i_Beg_stim = zeros(1,55);
i_End_stim = zeros(1,55);

is_real  = (intensities_order ~= 0);
cum_real = cumsum(is_real);
onsets   = triggers_onset(:).';
n_onsets = numel(onsets);
guard    = round(0.1*fs);    % 100 ms

for expo = 1:55
    inten = intensities_order(expo);
    if inten ~= 0
        k  = cum_real(expo);                 % k-th real
        k  = min(max(k,1), n_onsets);
        i0 = onsets(k);
    else
        k_prev = cum_real(expo);
        k_next = k_prev + 1;
        i0_after  = [];
        i0_before = [];
        if k_prev >= 1 && k_prev <= n_onsets
            i0_after  = onsets(k_prev) + round(sham_margin_sec*fs);
        end
        if k_next <= n_onsets
            i0_before = onsets(k_next) - round(sham_margin_sec*fs) - Lepoch + 1;
        end
        if expo == 1 || isempty(i0_after)
            i0 = i0_before;
        else
            i0 = i0_after;
            if ~isempty(i0_before)
                i1_after = i0_after + Lepoch - 1;
                guard_onset = onsets(min(k_next, n_onsets)) - guard;
                if i1_after >= guard_onset
                    i0 = i0_before;
                end
            end
        end
        if isempty(i0)
            i0 = max(1, N - Lepoch + 1);
        end
    end

    % Clamp to 1..N but DO NOT shift left to force full length.
    % (Short last epochs are handled later via padding.)
    i0 = max(1, i0);
    i1 = min(N, i0 + Lepoch - 1);

    i_Beg_stim(expo) = i0;
    i_End_stim(expo) = i1;
end
end