function [sep_2d, sep_3d] = calc_sep(array1, array2, pixel_size, z_step)
%CALC_SEP Calculate separation of coordinate arrays parsed from data_cell

%Convert the x, y, and z coords to nm
[xy1_nm, xyz1_nm] = nm_convert(array1, pixel_size, z_step);
[xy2_nm, xyz2_nm] = nm_convert(array2, pixel_size, z_step);
%Subtract and square the nm arays
xy_ss = (xy1_nm - xy2_nm).^2;
xyz_ss = (xyz1_nm - xyz2_nm).^2;
%Cacl separation in 2D and 3D
sep_2d = sqrt(xy_ss(:,1) + xy_ss(:,2));
sep_3d = sqrt(xyz_ss(:,1) + xyz_ss(:,2) + xyz_ss(:,3));
