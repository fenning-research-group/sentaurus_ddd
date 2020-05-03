# Load TDR file.
set mydata2D [load_file ${tdr_file}]

# Create new plot.
set myplot2D [create_plot -dataset $$mydata2D]

# Create 1D cutline normal to x-axis at point x=-0.005.
set mydata1D [create_cutline -plot $$myplot2D -type x -at ${x_cut_line}]

export_variables {${variables}} \
-dataset $$mydata1D -filename "${tdr_file}_band_diagram.csv" -overwrite