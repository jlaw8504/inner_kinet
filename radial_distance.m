function [rad_dist_2d, rad_dist_3d] = radial_distance(gfp1, gfp2, rfp1, rfp2, pixel_size, z_step)
%RADIAL_DISTANCE Returns the distance of each gfp signal from the spindle
%axis in both 2D and 3D

%Convert the x, y, and z coords to nm
[g1_xy, g1] = nm_convert(gfp1, pixel_size, z_step);
[g2_xy, g2] = nm_convert(gfp2, pixel_size, z_step);
[r1_xy, r1] = nm_convert(rfp1, pixel_size, z_step);
[r2_xy, r2] = nm_convert(rfp2, pixel_size, z_step);
%Calc 2D 3D radial distance
rad_dist_3d = [];
rad_dist_2d = [];
for n = 1:size(r1,1)
    rad_dist_3d(n,1) = norm(cross(r2(n,:)-r1(n,:),g1(n,:)-r1(n,:)))/...
        norm(r2(n,:)-r1(n,:));
    rad_dist_3d(n,2) = norm(cross(r2(n,:)-r1(n,:),g2(n,:)-r1(n,:)))/...
        norm(r2(n,:)-r1(n,:));
    rad_dist_2d(n,1) = abs(det([r2_xy(n,:)-r1_xy(n,:);g1_xy(n,:)-r1_xy(n,:)]))/...
        norm(r2_xy(n,:)-r1_xy(n,:));
    rad_dist_2d(n,2) = abs(det([r2_xy(n,:)-r1_xy(n,:);g2_xy(n,:)-r1_xy(n,:)]))/...
        norm(r2_xy(n,:)-r1_xy(n,:));    
end