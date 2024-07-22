# Use of diffTF with your own genome assembly and HOCOMOCO v12 TF Database

## DATA and tools used

- Genome assembly used: **mm39** found on ENSEMBL website: [Mus_musculus.GRCm39.dna_sm.toplevel](https://ftp.ensembl.org/pub/release-112/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz)

This choice was made base on the recomendation on bulk RNA-Seq pipeline from nf-core. The same genome assembly was used for the ATAC-Seq pipeline (nf-core).

- HOCOMOCOv12 logodds matrix. Database of TFs in mouse or Human , updated in 2024. Used by PWMScan

- PWMScan: to find TFs binding site on mm 39

- RNAseq results: counts matrix
- ATACSeq results: BAM file and peak annotation (MACS2)

- diffTF: it's a report on using diffTF, so obviously.

## Scan of TFBS based on HOCOMOCOv12

To use diffTF, you need the localisation of each potential binding site of your TFs of interest. They are bed file with the genomic coordinate, a score and the strand position (backward/forward). For more information, see the doc of [diffTF](https://difftf.readthedocs.io/en/latest/chapter2.html#input) on TFBS (transcription factor binding site). To create such file you will need the motif of each TF, i.e, a matrix of probability for the 4 nucleotides (MEME format) and a scanner. For that, I have used PWMScan

### Why using PWMScan?

1. You can use any genome assembly and TF database (here mm39 and HOCOMOCOv12)
2. Runs on your local computer (linux or mac os)
3. Run fast! ~10s/TFs (Linux with 12 core)
4. Easy to use
5. The output is exactly in the format needed for diffTF

In comparison, if you want to use FIMO you should wait more than 40 days for the 720 TFs of HOCOMOCOv12 (mouse) as it run on their server with a waiting list

### Prepare the Data for PWMScan

1. Clone the repo of PWM scan ([link](https://gitlab.sib.swiss/EPD/pwmscan)) and install it on your computer (linux or max os), see their documentation for installation

2. Download the logodds matrix of HOCOMOCO v12, [here](https://epd.expasy.org/ftp/pwmlib/hocomocov12_core_matrix_logodds.mat)

3. Create a unique file for each Matrix. You can do it with a simple bash script:

```bash
#!/bin/bash

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

```

#### NOTE:
The logodds matrix file will create 1443 unique file. There is human and mouse TFs. To filter the chosen one you can found all the informations needed on the HOCOMOCO website. The filtering of wanted TFs is let to your conveniance with any tool that you want (Python/R/bash/...). Also, for some TFs, there exist multiple possible motifs. In this report I use only the most probable motif so a total of 720 unique TFs for the mice.

### Use PWM Scan

Now that you have a logodds matrix for each TF you can use PWMScan. You will need two command: `matrix_prob` and `matrix_scan`. This two commands will be used for each TF individually as follow:
1. `matrix_prob -e 0.00001 -b 0.29,0.21,0.21,0.29 AHR.H12CORE.0.P.B.mat`
> return SCORE: 1305
2. `matrix_scan -m AHR.H12CORE.0.P.B.mat -c 1305 < mm39.fa > output_file.bed`

Here is a complete script for all TFs:

```bash
!/bin/bash
directory="data/matrix_example_to_scan"
genome="mm39.fa"

#check if the target is a directory
if [ ! -d "$directory" ]; then
	exit 1
fi


#loop through files in the directory
for file in "$directory"/*;
do
	echo "$file"
	output=$(./bin/matrix_prob -e 0.00001 -b 0.29,0.21,0.21,0.29 $file)
	score=$(echo "$output" | grep -oP 'SCORE\s*:\s*\K-?[0-9]+')
	echo "extracted SCORE: $score"
	file_name=$(basename "$file" .mat)
	./bin/matrix_scan -m $file -c $score <$genome > ../results/PWMScan_hocomocov12/$file_name.bed
done

echo "SCAN COMPLETED"
```

Here you should have 720 TFBS bed file (if you have previously filter your TF of interest).

#### NOTE
- The scan use a p-value threshold, here 0.00001. This can change a lot the number of possible binding site and so the output of diffTF. Check PWMScan documentation and diffTF documentation for more information.
- PWMScan output will directly be in the good format for diffTF so 3 columns for genomic coordinate, 1 column for the DNA motif, 1 column for the score and 1 column for the strand.
- You will need to rename the files, see point bellow

## Prepare the files in the good format for diffTF

Now you should have BED files that contains all the TFBS for each TFs of interest from the HOCOMOCOV12 database.

Before using diffTF you still need to refine your BED files in order to not crash the software. 

1. Rename BED file as followed: TFNAME_TFBS.bed
2. Filter unwanted sites
3. map your TF name to ENSEMBL name (for RNA-seq analysis)
4. Correct chromosome name. This point is only useful if your genome assembly doens't follow the convention to name their chromosome by using "chr" in front the chromosome number.

### 1. Rename BED files
to rename the bed files you can again use a bash script to rename all your files or any tools that you want.

```bash
#!/bin/bash

# Define the directory containing the .bed files
directory="../results/PWMScan_hocomocov12"

# Loop through all .bed files in the directory
for file in "$directory"/*.bed; do
    # Extract the base name
    base_name=$(basename ${file})
    # Use cut to get the NAME part
    name=$(echo "$base_name" | cut -d '.' -f 1)
    # Construct the new file name
    new_name="${directory}/${name}_TFBS.bed"
    # Rename the file
    mv "$file" "$new_name"
done
```

### 2. Filter unwanted genes
If your genome assembly has other location that standard chromosome like bellow, you will want to remove them in order to run correctly diffTF
```
JH584299.1	243370	243381	GTGCGTGCGTG	1278	+
JH584299.1	411755	411766	TTGCGTGCCTT	1285	-
JH584299.1	543109	543120	GTGCGTGCGTG	1278	+
JH584299.1	902889	902900	GTGCGTGCGTG	1278	+
GL456233.2	94368	94379	GTGCGTGCCTG	1317	-
JH584301.1	20638	20649	GTGCGTGCAAA	1313	+
GL456211.1	26774	26785	GTGCGTGCCAC	1360	+
GL456211.1	154506	154517	GTGCGTGCCAC	1360	+
GL456221.1	26692	26703	GTGCGTGCCAC	1360	+
GL456221.1	80466	80477	GTGCGTGCCAC	1360	-
GL456221.1	163896	163907	GTGCGTGCGTG	1278	+
JH584297.1	41270	41281	GTGCGTGCGTG	1278	-
JH584297.1	41274	41285	GTGCGTGCGTG	1278	-
JH584297.1	41278	41289	GTGCGTGCGTG	1278	-
JH584297.1	187815	187826	TTGCGTGCCTT	1285	+
JH584296.1	120886	120897	TTGCGTGCCTT	1285	+
JH584296.1	154503	154514	GTGCGTGCGTG	1278	-
GL456354.1	130378	130389	TTGCGTGCCTT	1285	+
JH584298.1	160636	160647	GTGCGTGCGTG	1278	+
JH584300.1	14213	14224	GTGCGTGCAAA	1313	+
JH584303.1	76391	76402	GTGCGTGCCAC	1360	-
GL456212.1	117894	117905	TTGCGTGAGAA	1368	+
GL456212.1	139860	139871	GTGCGTGCCAC	1360	+
MU069435.1	619	630	TCGCGTGCGAG	1307	-
MU069435.1	30747	30758	GTGCGTGCGTG	1278	+
```

you can remove them like that:

```bash
#!/bin/bash

# Define the directory containing the .bed files
directory="../results/PWMScan_hocomocov12"

# Loop through each .bed file in the directory
for file in "$directory"/*.bed; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Use grep to filter lines that begin with chr[1-19] or chr[X,Y]
        grep -E '^chr([1-9]|1[0-9]|X|Y)\b' "$file" > "${file}.tmp"
        
        # Overwrite the original file with the filtered content
        mv "${file}.tmp" "$file"
        
        echo "Processed $file"
    fi
done
```

### 3. Mapping TF name to ENSEMBL name

The file is a three column space separate file with SYMBOL (Name) ENSEMBL (ENSM0000X), HOCOID (NAME). I use the tool [g:Profiler](https://biit.cs.ut.ee/gprofiler/convert) to convert the name of my TFs and manual annotation for the few that were not found. There is probably better tool to do so, feel free to do it as you want.

Be carefull, the HOCOID are the name of your TFBS.bed file that can be completely different from the symbol name! So you need to have all of your TFBS.bed file with one entry in the HOCOID column.

#### TIPS
you can have a better name coverage with g:profiler by using the name from the json file from [HOCOMOCOv12](https://hocomoco12.autosome.org/downloads_v12) and not the TF name from the bed file.



### 4. Chromosome Name
This is probably the most important part of this report. diffTF work **uniquely** if your genome assembly fasta file use the convention to name the chromosome with "chr" at the beginning!! If, as in this example, your genome assembly begin with the chromosome number without the prefix "chr", you will need to change several files! 

Hopefully you will be able to run diffTF even if you have done your alignement with the "bad" convention. You will need to change several file for that

1. Add chr at the beginning of each line in all TFBS.bed file AND peak annotation
2. Add chr to your ATAC-Seq BAM file
3. Add chr to the FASTA file of your genome assembly


#### Add chr to TFBS.bed file

```bash
#!/bin/bash

# Define the directory containing the files
directory="PWMScan_HOCOMOCOv12"

# Loop through each file in the directory
for file in "$directory"/*; do
    # Check if the file is a regular file
    if [ -f "$file" ]; then
        # Use sed to add "chr" in front of every line and overwrite the file
        sed -i 's/^/chr/' "$file"
        echo "Updated $file"
    fi
done

```
if you have use MACS2 for your peak annotation, an output file .broadpeak or .narrowpeak will be created with a standard bed file for the beginning. this file also need a chr prefix in front of each line. You can use the same script as above.

#### Add chr to BAM file
This solution was found on biostars [here](https://www.biostars.org/p/13462/)
```bash
#!/bin/bash

for file in *.bam
do
  echo "$file"
  filename=`echo $file | cut -d "." -f 1`
  samtools view -H $file | \
      sed -e 's/SN:1/SN:chr1/' | sed -e 's/SN:2/SN:chr2/' | \
      sed -e 's/SN:3/SN:chr3/' | sed -e 's/SN:4/SN:chr4/' | \
      sed -e 's/SN:5/SN:chr5/' | sed -e 's/SN:6/SN:chr6/' | \
      sed -e 's/SN:7/SN:chr7/' | sed -e 's/SN:8/SN:chr8/' | \
      sed -e 's/SN:9/SN:chr9/' | sed -e 's/SN:10/SN:chr10/' | \
      sed -e 's/SN:11/SN:chr11/' | sed -e 's/SN:12/SN:chr12/' | \
      sed -e 's/SN:13/SN:chr13/' | sed -e 's/SN:14/SN:chr14/' | \
      sed -e 's/SN:15/SN:chr15/' | sed -e 's/SN:16/SN:chr16/' | \
      sed -e 's/SN:17/SN:chr17/' | sed -e 's/SN:18/SN:chr18/' | \
      sed -e 's/SN:19/SN:chr19/' | sed -e 's/SN:20/SN:chr20/' | \
      sed -e 's/SN:21/SN:chr21/' | sed -e 's/SN:22/SN:chr22/' | \
      sed -e 's/SN:X/SN:chrX/' | sed -e 's/SN:Y/SN:chrY/' | \
      sed -e 's/SN:MT/SN:chrM/' | samtools reheader - $file > ${filename}_chr.bam
done

```

#### Add chr to FASTA file

```bash
sed -f pattern.txt < Mus_musculus.GRCm39.dna_sm.primary_assembly.fa > mm39.fa

```
with the file pattern.txt:

s/^>1/>chr1/
s/^>2/>chr2/
s/^>3/>chr3/
s/^>4/>chr4/
s/^>5/>chr5/
s/^>6/>chr6/
s/^>7/>chr7/
s/^>8/>chr8/
s/^>9/>chr9/
s/^>10/>chr10/
s/^>11/>chr11/
s/^>12/>chr12/
s/^>13/>chr13/
s/^>14/>chr14/
s/^>15/>chr15/
s/^>16/>chr16/
s/^>17/>chr17/
s/^>18/>chr18/
s/^>19/>chr19/
s/^>X/>chrX/
s/^>Y/>chrY/

### Run the analysis
Running diff TF is very easy if you use singularity (as recommanded by the authors)! The installation of diffTF and utilisation will not be cover here as the documentation is already well detailed, See [documentation](https://difftf.readthedocs.io/en/latest/index.html)

The thing that you will need to update with your own genome assembly and new TFBS in addition to your sample data (see documentation) are listed here:

- In the confiFile.json
    1. refGenome_fasta: path/to/your/genomeAssembly.fa
    2. dir_TFBS: path/to/your/TFBS.bed 
    3. HOCOMOCO_mapping: path/to/your/translationTable.csv (the SYMBOL/ENSM/HOCOID file)
1. Use PWMScan, see github
2. Change name to NAME_TFBS.bed
3. Add chr in front of every files
4. create a spreadsheet to HOCOMOCO name and ENS ID and gene Symbol --> use https://biit.cs.ut.ee/gprofiler/convert
