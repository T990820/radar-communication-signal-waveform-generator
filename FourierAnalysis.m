function [] = FourierAnalysis(path)
s = load(path);
s = s.('s');
F = fft(s);
F = abs(F) / max(abs(F));
s = F;
save(path,"s");
end

