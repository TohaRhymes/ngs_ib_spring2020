#illumina assembling

spades -1 data/illumina.100x.1.fq -2 data/illumina.100x.2.fq -o "illumina_assembling"



# Hybrid assembling

spades --pacbio data/pacbio_10x.fq -1 data/illumina.100x.1.fq -2 data/illumina.100x.2.fq -o "hybrid_assembling_cov10"

spades --pacbio data/pacbio_20x.fq -1 data/illumina.100x.1.fq -2 data/illumina.100x.2.fq -o "hybrid_assembling_cov20"

spades --pacbio data/pacbio_40x.fq -1 data/illumina.100x.1.fq -2 data/illumina.100x.2.fq -o "hybrid_assembling_cov40"

spades --pacbio data/pacbio_80x.fq -1 data/illumina.100x.1.fq -2 data/illumina.100x.2.fq -o "hybrid_assembling_cov80"
