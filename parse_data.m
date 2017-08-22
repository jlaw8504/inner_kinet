function [gfp1, gfp2, rfp1, rfp2, pixel_size, z_step] = parse_data(directory)
%PARSE_DATA Collect coordinate data from heatmap_GUI.m .mat output files
%Instantiate arrays
gfp1 = [];
gfp2 = [];
rfp1 = [];
rfp2 = [];
pixel_size = [];
z_step = [];
%Loop through all .mat files in provided directory
files = dir(fullfile(directory,'*.mat'));
for n = 1:size(files,1)
    load(fullfile(directory,files(n).name));
    %parse the data_cell cell structure
    for i = 2:size(data_cell,1)
        gfp1 = [gfp1; data_cell{i,1}];
        gfp2 = [gfp2; data_cell{i,3}];
        rfp1 = [rfp1; data_cell{i,5}];
        rfp2 = [rfp2; data_cell{i,6}];
        pixel_size = [pixel_size; data_cell{i,8}];
        z_step = [z_step; data_cell{i,7}];
    end
end