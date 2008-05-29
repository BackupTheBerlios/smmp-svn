bg_color white
load best.pdb
orient
show cartoon
set ray_trace_mode, 0
set ray_trace_fog = 0
set antialias, 2
ray 320, 240
png best.png