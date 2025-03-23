function xy = mni2xy_v2(x, y, z, view)

x_scale = 0.87;
y_scale = 0.95;

switch view
    case 'sagittal'
        xy(:,1) = x_scale*y+98;
        xy(:,2) = y_scale*z+110;
    case 'axial'
        xy(:,1) = x_scale*x + 8;
        xy(:,2) = y_scale*y + 2;
    case 'coronal'
        xy(:,1) = x_scale*x + 8;
        xy(:,2) = y_scale*z+112;
end