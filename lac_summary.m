function s = lac_summary(directory, slbounds)
%LAC_SUMMARY Summary script for comparing WT and brn1-9 mutants containing
%1p7 kb lacO/LacI-GFP array

%Parse data
[s.gfp1, s.gfp2, s.rfp1, s.rfp2, s.pixel_size, s.z_step] = parse_data(directory);
%Spindle Length
[s.sl_2d, s.sl_3d] = calc_sep(s.rfp1, s.rfp2, s.pixel_size, s.z_step);

%Lac Separation
[s.sep_2d, s.sep_3d] = calc_sep(s.gfp1, s.gfp2, s.pixel_size, s.z_step);
s.mean_sep_2d = mean(s.sep_2d(:));
s.sem_sep_2d = std(s.sep_2d(:))/sqrt(numel(s.sep_2d(:)));
s.mean_sep_3d = mean(s.sep_3d(:));
s.sem_sep_3d = std(s.sep_3d(:))/sqrt(numel(s.sep_3d(:)));
%Radial displacement analysis
[s.rad_dist_2d, s.rad_dist_3d] = radial_distance(s.gfp1, s.gfp2, s.rfp1, s.rfp2, s.pixel_size, s.z_step);
s.mean_rd2 = mean(s.rad_dist_2d(:));
s.sem_rd2 = std(s.rad_dist_2d(:))/sqrt(numel(s.rad_dist_2d(:)));
s.mean_rd3 = mean(s.rad_dist_3d(:));
s.sem_rd3 = std(s.rad_dist_3d(:))/sqrt(numel(s.rad_dist_3d(:)));