function [xy, xyz] = nm_convert(array, pixel_size, z_step)
%NM_CONVERT Convert array [x, y, z, intensity] parsed from data cell into
%two arrays, [x, y] and [z] in nm.

%Convert the x and y coords to nm using pixel_size array
xy = array(:,1:2) .* [pixel_size,pixel_size];
%Convert the plane to nm using z_step
z = array(:,3) .* z_step;
%Concatenate xy and z for three dimensional matrix
xyz = [xy,z];
