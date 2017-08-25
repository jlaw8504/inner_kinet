function [gfp1, gfp2, rfp1, rfp2, pixel_size, z_step, stretch_array] = spindle_bounds(directory, slbounds, dim)
%SPINDLE_BOUNDS Calculates the spindle length of each cell in 2D or 3D.
%Restricts the output arrays to only data that falls within spindle bounds

%Run parse_data to gather initial arrays, pixel size, and z-step size
[gfp1, gfp2, rfp1, rfp2, pixel_size, z_step, stretch_array] = parse_data(directory);

%Run calc_sep function to determine spindle lengths in 2D and 3D
[sl_2d, sl_3d] = calc_sep(rfp1, rfp2, pixel_size, z_step);
%Create binary arrays based on spindle-length bounds
if dim == 2
    sl_bin = sl_2d > slbounds(1) & sl_2d < slbounds(2);
elseif dim == 3
    sl_bin = sl_3d > slbounds(1) & sl_3d < slbounds(2);
else
    error('Please use 2 or 3 as dim variable');
end

%Filter coordinate arrays pixel_size, and z-step using sl_bin
gfp1 = gfp1(sl_bin,:);
gfp2 = gfp2(sl_bin,:);
rfp1 = rfp1(sl_bin,:);
rfp2 = rfp2(sl_bin,:);
pixel_size = pixel_size(sl_bin,:);
z_step = z_step(sl_bin,:);
stretch_array = stretch_array(sl_bin,:);

