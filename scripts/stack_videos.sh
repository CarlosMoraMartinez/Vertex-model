
ffmpeg -i hexgrid2021_19_0_21grid5.avi -i hexgrid2021_19_1_21grid5.avi -i hexgrid2021_19_2_21grid5.avi -i hexgrid2021_19_3_21grid5.avi -i hexgrid2021_19_4_21grid5.avi \
-filter_complex \
"[0:v][1:v]hstack=inputs=2[top];\
[2:v][3:v]hstack=inputs=2[bottom];\
[top][bottom]vstack=inputs=2[v]" -map "[v]" output.avi



ffmpeg -i hexgrid2021_19_0_21grid5.avi -i hexgrid2021_19_1_21grid5.avi -i hexgrid2021_19_2_21grid5.avi -i hexgrid2021_19_3_21grid5.avi -i hexgrid2021_19_4_21grid5.avi -i hexgrid2021_19_5_21grid5.avi -i hexgrid2021_19_6_21grid5.avi -i hexgrid2021_19_7_21grid5.avi -i hexgrid2021_19_8_21grid5.avi \
-filter_complex \
"[0:v][1:v][2:v]hstack=inputs=3[top];\
[3:v][4:v][5:v]hstack=inputs=3[middle];\
[6:v][7:v][8:v]hstack=inputs=3[bottom];\
[top][middle][bottom]vstack=inputs=3[v]" -map "[v]" output2.avi




ffmpeg -i hexgrid2021_17_6_21grid5.avi -i hexgrid2021_17_7_21grid5.avi -i hexgrid2021_17_8_21grid5.avi -i hexgrid2021_17_9_21grid5.avi -i hexgrid2021_17_10_21grid5.avi -i hexgrid2021_17_11_21grid5.avi -i hexgrid2021_17_12_21grid5.avi -i hexgrid2021_17_13_21grid5.avi -i hexgrid2021_17_14_21grid5.avi \
-filter_complex \
"[0:v][1:v][2:v]hstack=inputs=3[top];\
[3:v][4:v][5:v]hstack=inputs=3[middle];\
[6:v][7:v][8:v]hstack=inputs=3[bottom];\
[top][middle][bottom]vstack=inputs=3[v]" -map "[v]" output2.avi
