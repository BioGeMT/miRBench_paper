To create a text file with a list of fastq filenames in a directory (excluding files ending with '_2.fastq.gz'), run:
`bash get_filenames.sh <source directory> <output file>`

To create symlinks in a destination directory, to files in a source directory (excluding files ending with '_2.fastq.gz'), run:
`bash create_symlinks.sh <source directory> <destination directory>`