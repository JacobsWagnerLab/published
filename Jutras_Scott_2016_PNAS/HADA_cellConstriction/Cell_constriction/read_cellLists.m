%Core script for loading experimental data

%location of cellLists
load('cp32_GFP_HADA_DRAQ5_meshes.mat')
disp('loading images')
%corresponding locations of images
dna = loadIMstack('cp32_GFP_HADA_DRAQ5_dna.tif');
gfp = loadIMstack('cp32_GFP_HADA_DRAQ5_gfp.tif');
hada = loadIMstack('cp32_GFP_HADA_DRAQ5_hada.tif');
phs = loadIMstack('cp32_GFP_HADA_DRAQ5_phase.tif');
disp('checking for saturated cells and calculating signals')
ims = {dna,gfp,hada};
[cellList,~] = removeSaturatedCells(cellList, ims);
signals = reformatCellList(cellList);

disp('loading cellList')
load('cp32_GFP_HADA_DRAQ5_0331_meshes.mat')
disp('loading images')
phase = loadIMstack('cp32_GFP_HADA_DRAQ5_0331_phase.tif');
gfp = loadIMstack('cp32_GFP_HADA_DRAQ5_0331_gfp.tif');
dna = loadIMstack('cp32_GFP_HADA_DRAQ5_0331_dna.tif');
hada = loadIMstack('cp32_GFP_HADA_DRAQ5_0331_hada.tif');
disp('checking for saturated cells and calculating signals')
ims = {dna,gfp,hada};
[cellList,~] = removeSaturatedCells(cellList, ims);

signals = [signals, reformatCellList(cellList)];

save('cellList_signals','signals')