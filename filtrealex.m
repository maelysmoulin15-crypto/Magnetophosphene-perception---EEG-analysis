function res=filtrealex (fname, sr, hi, ramp1, lo, ramp2)
% sr = Sampling rate (ex: 100)
% hi = highpass frequency cut (ex: 3Hz)
% lo = lowpass frequency cut (ex: 20Hz)
% Filtre passe haut, passe bas


fname = fname - mean(fname);

hp = 1-makefilt(hi/sr,length(fname),ramp1/sr);

hfname = datafilt (fname, hp);

low1520 = makefilt (lo/sr, length(hfname), ramp2/sr);

filtdata = datafilt (hfname,low1520);

% assignin ('base','filtdata',filtdata)

res = filtdata;


return
