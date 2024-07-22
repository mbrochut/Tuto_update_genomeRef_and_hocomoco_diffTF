# Define the input file and output directory
input_file="../data/hocomocov12_core_matrix_logodds.mat"
output_dir="../results/matrix_hocomoco_v12_split"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Use awk to split the file into separate .mat files
awk -v output_dir="$output_dir" '
BEGIN { file_count = 0; }
/^>log-odds matrix/ {
    if (file_count > 0) {
        close(output_file);
    }
    file_count++;
    split($0, name, " ");
    # Remove colons from the name[3] part
    gsub(":", "", name[3]);
    output_file = output_dir "/" name[3] ".mat";
}
{
    print > output_file;
}
' "$input_file"

echo "Splitting completed."