This script filters .tsv files from HybriDetector for miRNA sequences, 
extracts specific information (id, seq.g, seq.m, noncodingRNA_fam, feature,test,label) about the miRNA sequences only 
and creates a new .tsv file containing the above information

To use:

"python miRNAfinder.py input_file output_file"

If the output file is not specified, the output is saved in a file named "filtered_data.tsv" in the current directory. If more 2 arguments are given, they are ignored.


EXAMPLE INPUT:

chr.g	start.g	end.g	flag.g	strand.g	gene_name	gene_biotype	feature	seq.g	noncodingRNA_real_seq	noncodingRNA	noncodingRNA_type	noncodingRNA_fam	forward_dir	Ndups	Nunique	dups_perc	repeatmasker	repeat_fam	noncodingRNA_nomm	mrna_length	cov	Nreads	mrna_aln	mrna_aln2	mrna_aln5	mrna_aln10	mrna_aln20	A_perc	C_perc	G_perc	T_perc	mfe_struct	mfe	ene	noncodingRNA_seq
1	7779809	7779858	0	+	VAMP3	protein_coding	three_prime_utr	TCCAAATTAGAAGGCCGGCCCCGTCCACATTTTGCACAGTGCCTTTACAG	TGTGCAAATCCATGCAAAACTGA	hsa-miR-19b-3p|hsa-miR-19a-3p	miRNA	mir-19	FALSE	101	56	44.6			TRUE	TRUE	785.54	5475	TRUE	TRUE	TRUE	TRUE	TRUE	39.1	21.7	17.4	21.7	............((.((....)).))....((((((((............&))))))))...............	-12.8	1.8	TGTGCAAATCCATGCAAAACTGA
8	37698821	37698870	0	+	ZNF703	protein_coding	three_prime_utr	TATTTACTATAATGTTAGCTTACAAGCTGGGAATATAAGTGCATTAACGG	TATTGCACTTGTCCCGGCCTGT	hsa-miR-92a-3p|hsa-miR-92b-3p	miRNA	mir-25	FALSE	39	28	28.2			TRUE	TRUE	2106.973	14685	TRUE	TRUE	TRUE	TRUE	TRUE	9.1	31.8	22.7	36.4	......(((......)))...(((.(((((((...((((((((.((....&)).))))))))))))))).)))	-24.9	5	TATTGCACTTGTCCCGGCCTGT
7	137879851	137879900	16	-	CREB3L2	protein_coding	three_prime_utr	ATATTCTCAGGCCGCAAGTGCAATGCCTGAGGGGATCAGGCTTTTCTACT	TATTGCACTTGTCCCGGCCTGT	hsa-miR-92a-3p|hsa-miR-92b-3p	miRNA	mir-25	FALSE	32	21	34.4			TRUE	TRUE	323.5425	2255	TRUE	TRUE	TRUE	TRUE	TRUE	9.1	31.8	22.7	36.4	.......(((((((((((((((((((((((.....)))))).........&.))))))))))...))))))).	-37.2	5	TATTGCACTTGTCCCGGCCTGT
10	58767584	58767633	16	-	BICC1	protein_coding	intron	AAAGCAATTATTTTTAAATTCCTTGTCGGGTAATTTGTAGACCTCTACTT	GGTCCTAAGGTAGCGA	12_127619336_127619447_-_Y_RNA.49-201_10	YRNA	0	FALSE	29	18	37.9	L1M3c	LINE	TRUE	TRUE	5.59563	39	TRUE	TRUE	TRUE	TRUE	TRUE	25	18.8	37.5	18.8	..............(((((((((....))).))))))..((((.......&))))............	-9.4	0	GGCTGGTCCTAAGGTAGCGAGTTATCTCAATCGACTGTTCACAGTCAGTTACAGATCAACCTCCTTGTTCTACTCTTTCCCTTCTTCCACTACTCCACTTGACTAGTCTAAA
16	3022076	3022125	0	+	TNFRSF12A	protein_coding	three_prime_utr	ACTAAGGAACTGCAGCATTTGCACAGGGGAGGGGGGTGCCCTCCTTCCTA	TGTGCAAATCCATGCAAAACTGA	hsa-miR-19b-3p|hsa-miR-19a-3p	miRNA	mir-19	FALSE	30	18	40			TRUE	TRUE	366.5858	2555	TRUE	TRUE	TRUE	TRUE	TRUE	39.1	21.7	17.4	21.7	..........((((..(((((((((((((((((((....)))))))))).&)))))))))...)))).......	-33.1	5	TGTGCAAATCCATGCAAAACTGA
X	151408990	151409039	0	+	VMA21	protein_coding	three_prime_utr	TGCACAGCACCTTACAGTTTGCAAAGAACGTTCACACATTCTCATTTGAG	TAAGGTGCATCTAGTGCAGATAG	hsa-miR-18a-5p	miRNA	mir-17	FALSE	27	17	37	MIRc	SINE	TRUE	TRUE	994.301	6930	TRUE	TRUE	TRUE	TRUE	TRUE	30.4	13	30.4	26.1	(((((.((((((((..((((.....))))...........(((....)))&)))))))).....))))).....	-18.4	5	TAAGGTGCATCTAGTGCAGATAG

EXAMPLE OUTPUT:

id	seq.g	seq.m	noncodingRNA_fam	feature	test	label
1	TCCAAATTAGAAGGCCGGCCCCGTCCACATTTTGCACAGTGCCTTTACAG	TGTGCAAATCCATGCAAAACTGA	mir-19	three_prime_utr	True	1
2	TATTTACTATAATGTTAGCTTACAAGCTGGGAATATAAGTGCATTAACGG	TATTGCACTTGTCCCGGCCTGT	mir-25	three_prime_utr	False	1
3	ATATTCTCAGGCCGCAAGTGCAATGCCTGAGGGGATCAGGCTTTTCTACT	TATTGCACTTGTCCCGGCCTGT	mir-25	three_prime_utr	False	1
4	ACTAAGGAACTGCAGCATTTGCACAGGGGAGGGGGGTGCCCTCCTTCCTA	TGTGCAAATCCATGCAAAACTGA	mir-19	three_prime_utr	False	1
5	TGCACAGCACCTTACAGTTTGCAAAGAACGTTCACACATTCTCATTTGAG	TAAGGTGCATCTAGTGCAGATAG	mir-17	three_prime_utr	False	1
6	TGTTTTGCACGATGTAAAAACCCTGTCTTTTTGCACGATACAGCCAAAAG	TGTGCAAATCCATGCAAAACTGA	mir-19	three_prime_utr	True	1
