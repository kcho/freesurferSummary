save_tiff $env(_SNAPSHOT_STEM)-lat.tif
make_lateral_view
rotate_brain_y 180
redraw
save_tiff $env(_SNAPSHOT_STEM)-med.tif
make_lateral_view
rotate_brain_x 90
redraw
save_tiff $env(_SNAPSHOT_STEM)-ven.tif
make_lateral_view
rotate_brain_x -90
redraw
save_tiff $env(_SNAPSHOT_STEM)-dor.tif
exit