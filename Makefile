DIR_DATA=../data07
FILE_REF=../ref/true_genes_only.gtf
FILE_REF_ALL=../ref/genomic.gtf
FILE_FNA=../ref/GCF_000750555.1_ASM75055v1_genomic.fna
DIR_COUNT_07="../data07/deseq2"
DIR_COUNT_08="../data08B/deseq2"

rDiff_PATH=/Users/deng/Research/rna/RiboDiff/scripts/TE.py

DIR_BOWTIE=$(DIR_DATA)/bowtie
DIR_COUNT=$(DIR_DATA)/feature
DIR_DESEQ2=$(DIR_DATA)/deseq2
DIR_VENN=$(DIR_DATA)/venn
DIR_NOISE=$(DIR_DATA)/noise
DIR_RNA=$(DIR_DATA)/feature_RNA
DIR_RDIFF=$(DIR_DATA)/rdiff
DIR_RIBO_VENN=$(DIR_DATA)/venn_ribo
DIR_ASD=$(DIR_DATA)/ASD

DIR_MOTIF=$(DIR_DATA)/motif

all: 

count:$(DIR_COUNT)

$(DIR_COUNT): 2-count.sh $(DIR_BOWTIE)
	mkdir -p $(DIR_COUNT)
	./2-count.sh $(FILE_REF) $(DIR_BOWTIE) $(DIR_COUNT)

deseq2:$(DIR_DESEQ2)

$(DIR_DESEQ2):  3-to_csv_feature.py 4-deseq2.R
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


RNA: $(DIR_RNA)


$(DIR_RNA): ribo-count.sh $(FILE_REF_rRNA) $(DIR_BOWTIE)
	mkdir -p $(DIR_RNA)
	./ribo-count.sh $(FILE_REF_ALL) $(DIR_BOWTIE) $(DIR_RNA)


rDiff:
	mkdir -p $(DIR_RDIFF) 
	python ribo_input.py 0 $(DIR_COUNT_07) $(DIR_COUNT_08) $(DIR_RDIFF)
	#source init.sh
	#conda init bash
	#conda deactivate
	#conda activate ribo
	/Users/deng/.conda/envs/ribo/bin/python $(rDiff_PATH) -e $(DIR_RDIFF)/meta-1.csv -c $(DIR_RDIFF)/matrix-1.csv -o $(DIR_RDIFF)/result-1.tsv -p 1
	#conda deactivate

rDiff-C:
	mkdir -p $(DIR_RDIFF) 
	python ribo_input.py 1 $(DIR_COUNT_07) $(DIR_COUNT_08) $(DIR_RDIFF)
	#init.sh
	/Users/deng/.conda/envs/ribo/bin/python $(rDiff_PATH) -e $(DIR_RDIFF)/meta-2.csv -c $(DIR_RDIFF)/matrix-2.csv -o $(DIR_RDIFF)/result-2.tsv -p 1

ribo-venn:
	mkdir -p $(DIR_RIBO_VENN)
	python venn.py $(DIR_RDIFF) $(DIR_RIBO_VENN)




########################

all-motif: seqs motif html

seqs: 
	mkdir -p $(DIR_MOTIF)
	python A1-meme-seq.py "$(DIR_RIBO_VENN)/B_m &C_p.txt" "$(DIR_MOTIF)/B_mC_p.fna" $(FILE_FNA) $(FILE_REF_ALL)
	python A1-meme-seq.py "$(DIR_RIBO_VENN)/B_p & C_m.txt" "$(DIR_MOTIF)/B_pC_m.fna" $(FILE_FNA) $(FILE_REF_ALL)
	python A1-meme-seq.py "$(DIR_RIBO_VENN)/B_m.txt" "$(DIR_MOTIF)/B_m.fna" $(FILE_FNA) $(FILE_REF_ALL)

motif:
	#meme "$(DIR_MOTIF)/B_mC_p.fna" -oc "$(DIR_MOTIF)/B_mC_p" -p 4 -dna -mod zoops -nmotifs 5
	#meme "$(DIR_MOTIF)/B_pC_m.fna" -oc "$(DIR_MOTIF)/B_pC_m" -p 4 -dna -mod zoops -nmotifs 5
	meme "$(DIR_MOTIF)/B_m.fna" -oc "$(DIR_MOTIF)/B_m" -p 4 -dna -mod zoops -nmotifs 5


# Prepare B_m with backgroud reference (<1.5 fold)
Bmbg:
	python venn_bg1.py $(DIR_RDIFF) $(DIR_RIBO_VENN)
	python A1-meme-seq.py "$(DIR_RIBO_VENN)/bg4Bm.txt" "$(DIR_MOTIF)/bg4Bm.fna" $(FILE_FNA) $(FILE_REF_ALL)
	meme "$(DIR_MOTIF)/B_m.fna" -objfun de -neg "$(DIR_MOTIF)/bg4Bm.fna" -oc "$(DIR_MOTIF)/Bm_with_bg" -p 4 -dna -mod zoops -nmotifs 5

# Prepare B_m with backgroud reference (<1 fold)
Bmbg2:
	python venn_bg.py $(DIR_RDIFF) $(DIR_RIBO_VENN)
	python A1-meme-seq.py "$(DIR_RIBO_VENN)/bg4Bm2.txt" "$(DIR_MOTIF)/bg4Bm2.fna" $(FILE_FNA) $(FILE_REF_ALL)
	meme "$(DIR_MOTIF)/B_m.fna" -objfun de -neg "$(DIR_MOTIF)/bg4Bm2.fna" -oc "$(DIR_MOTIF)/Bm_with_bg2" -p 4 -dna -mod zoops -nmotifs 5


TIR:
	python A2-TIR.py $(DIR_RIBO_VENN) "$(DIR_DATA)/TIR" $(FILE_FNA) $(FILE_REF_ALL)



TIR-fisher:
	python A2-fisher.py $(DIR_RIBO_VENN) "$(DIR_DATA)/TIR" $(FILE_FNA) $(FILE_REF_ALL) #> $(DIR_DATA)/TIR/t-test.txt

ASD:
	mkdir -p $(DIR_ASD)
	python C1-SD.py $(DIR_RIBO_VENN) $(DIR_ASD) $(FILE_FNA) $(FILE_REF_ALL) 


ASD-only:
	mkdir -p $(DIR_ASD)
	python C2-SD-only.py $(DIR_RIBO_VENN) $(DIR_ASD) $(FILE_FNA) $(FILE_REF_ALL) 

DIR_GO=$(DIR_DATA)/GO
GO:
	mkdir -p $(DIR_GO)
	python D1-go.py $(DIR_RIBO_VENN) $(DIR_GO) $(FILE_REF_ALL) 