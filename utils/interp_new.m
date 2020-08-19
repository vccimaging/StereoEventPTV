function interpimg = interp_new( x1, y1, z1, x ,x1_tmp, y1_tmp, z1_tmp, interp_mode)
interpimg =  interp3( x1, y1, z1, x, x1_tmp, y1_tmp, z1_tmp, interp_mode);
interpimg(isnan(interpimg)) = 0;
end