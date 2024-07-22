!/bin/bash
directory="../data/matrix_example_to_scan"
genome="../data/mm39.fa"

#check if the target is a directory
if [ ! -d "$directory" ]; then
	exit 1
fi


#loop through files in the directory
for file in "$directory"/*;
do
	echo "$file"
	output=$(./../pwmscan/bin/matrix_prob -e 0.00001 -b 0.29,0.21,0.21,0.29 $file)
	score=$(echo "$output" | grep -oP 'SCORE\s*:\s*\K-?[0-9]+')
	echo "extracted SCORE: $score"
	file_name=$(basename "$file" .mat)
	./../pwmscan/bin/matrix_scan -m $file -c $score <$genome > ../results/PWMScan_hocomocov12/$file_name.bed
done

echo "SCAN COMPLETED"