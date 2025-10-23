% === Auto-sync figures to Overleaf ===
src = 'Figures/';
dst = '/Users/sakiclaudia/Library/CloudStorage/Dropbox/Apps/Overleaf/JET_Portfolio_ConvenienceYield/NewCode/Figures/';

if ~exist(dst, 'dir')
    mkdir(dst);
end

extensions = {'*.pdf', '*.eps'};

for k = 1:length(extensions)
    files = dir(fullfile(src, extensions{k}));
    for i = 1:length(files)
        copyfile(fullfile(src, files(i).name), fullfile(dst, files(i).name));
    end
end

disp('Figures copied to Overleaf.');