function idx = pick_detection_channel(eeg)
% pick_detection_channel
% Selects the cleanest possible detection channel.
% Priority: M1/TP9/A1/LM/MASTOIDL… -> else MF/Laser -> else 65 -> else last.
labels = upper(string({eeg.chanlocs.labels}));
prefer = ["M1","TP9","A1","LM","LMASTOID","MASTOIDL","M1-","-M1","MASTOID L","LEFT MASTOID","A1-"];
for p = prefer
    hit = find(contains(labels, p), 1);
    if ~isempty(hit), idx = hit; return; end
end
mf_keys = ["LASER","MF","MAG","COIL","STIM","EXPOSURE"];
hit = find(ismember(labels, mf_keys), 1);
if ~isempty(hit), idx = hit; return; end
nCh = size(eeg.data,1);
if nCh >= 65, idx = 65; else, idx = nCh; end
end
