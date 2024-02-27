function [] = FourierAnalysis(folder)
paths = getAllFiles(folder);
for i = 1:size(paths,1)
    path = paths{i};
    s = load(path);
    s = s.('s');
    F = fft(s);
    F = abs(F) / max(abs(F));
    s = F;
    save(path,"s");
end
end

