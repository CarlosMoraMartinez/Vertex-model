ffmpeg -framerate 3 -i displacement4_displacementVectors_%d.png -codec copy displ4.avi
ffmpeg -framerate 3 -i displacement3_displacementVectors_%d.png -codec copy displ3.avi
ffmpeg -framerate 3 -i displacement2_displacementVectors_%d.png -codec copy displ2.avi
ffmpeg -framerate 3 -i displacement1_displacementVectors_%d.png -codec copy displ1.avi

