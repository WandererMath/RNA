DIR_DATA=../data08
FILE_REF=../ref/genes_only.gtf
FILE_REF_rRNA=../ref/rRNA-only.gtf
DIR_COUNT_07="../data07/deseq2"
DIR_COUNT_08="../data08/deseq2"

DIR_BOWTIE=$(DIR_DATA)/bowtie
DIR_COUNT=$(DIR_DATA)/feature
DIR_DESEQ2=$(DIR_DATA)/deseq2
DIR_VENN=$(DIR_DATA)/venn
DIR_NOISE=$(DIR_DATA)/noise
DIR_rRNA=$(DIR_DATA)/rRNA
DIR_RDIFF=$(DIR_DATA)/rdiff

all: $(DIR_DESEQ2)



$(DIR_COUNT): 2-count.sh $(DIR_BOWTIE)
	mkdir -p $(DIR_COUNT)
	./2-count.sh $(FILE_REF) $(DIR_BOWTIE) $(DIR_COUNT)

$(DIR_DESEQ2): $(DIR_COUNT) 3-to_csv_feature.py 4-deseq2.R
	mkdir -p $(DIR_DESEQ2)
	python 3-to_csv_feature.py 0 $(DIR_COUNT) $(DIR_DESEQ2)
	python 3-to_csv_feature.py 1 $(DIR_COUNT) $(DIR_DESEQ2)

	Rscript 4-deseq2.R $(DIR_DESEQ2) 0
	mv $(DIR_DESEQ2)/Rplots.pdf $(DIR_DESEQ2)/A-B.pdf
	Rscript 4-deseq2.R $(DIR_DESEQ2) 1
	mv $(DIR_DESEQ2)/Rplots.pdf $(DIR_DESEQ2)/A-C.pdf

$(DIR_DATA)/scatter/norm.csv: $(DIR_DATA)/feature 15-csv.py 16-norm.R 
	mkdir -p $(DIR_DATA)/scatter
	python 15-csv.py $(DIR_DATA)/feature $(DIR_DATA)/scatter
	Rscript 16-norm.R $(DIR_DATA)/scatter

scatter: $(DIR_DATA)/scatter/norm.csv
	python 17-corr.py $(DIR_DATA)/scatter

venn: $(DIR_DESEQ2) 14-venn_padj.py
	mkdir -p $(DIR_VENN)
	python 14-venn_padj.py $(DIR_DESEQ2) $(DIR_VENN)

noise:  17-noise.sh $(DIR_BOWTIE)
	mkdir -p $(DIR_NOISE)
	./17-noise.sh $(FILE_REF) $(DIR_BOWTIE) $(DIR_NOISE)

rRNA: $(DIR_rRNA)

$(DIR_rRNA): ribo-count.sh $(FILE_REF_rRNA) $(DIR_BOWTIE)
	mkdir -p $(DIR_rRNA)
	./ribo-count.sh $(FILE_REF_rRNA) $(DIR_BOWTIE) $(DIR_rRNA)

rDiff:
	mkdir -p $(DIR_RDIFF) 
	python ribo_input.py 0 $(DIR_COUNT_07) $(DIR_COUNT_08) $(DIR_RDIFF)
	#source init.sh
	#conda init bash
	#conda deactivate
	#conda activate ribo
	/Users/deng/.conda/envs/ribo/bin/python $(rDiff_PATH) -e $(DIR_RDIFF)/meta-1.csv -c $(DIR_RDIFF)/matrix-1.csv -o $(DIR_RDIFF)/result-1.tsv
	#conda deactivate

rDiff-C:
	mkdir -p $(DIR_RDIFF) 
	python ribo_input.py 1 $(DIR_COUNT_07) $(DIR_COUNT_08) $(DIR_RDIFF)

	/Users/deng/.conda/envs/ribo/bin/python $(rDiff_PATH) -e $(DIR_RDIFF)/meta-2.csv -c $(DIR_RDIFF)/matrix-2.csv -o $(DIR_RDIFF)/result-2.tsv


