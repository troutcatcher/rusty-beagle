#!/bin/bash
# Fetches the public synbreedData cattle dataset (Wimmer et al., GPL-2):
# 500 real dairy bulls genotyped at 7,250 SNPs across 29 autosomes, with
# pedigree. Source: CRAN package synbreedData 1.5 via the CRAN GitHub mirror.
# Requires R (r-base-core) to decode the RData file into TSVs.
set -eu
DIR=${1:-.}
URL=https://raw.githubusercontent.com/cran/synbreedData/master/data/cattle.RData
MD5=5eab41d298545576198487501fd99cad
cd "$DIR"
[ -f cattle.RData ] || curl -sSL -o cattle.RData "$URL"
echo "$MD5  cattle.RData" | md5sum -c -
Rscript - <<'RS'
load("cattle.RData")
g <- cattle$geno
num <- matrix(NA_integer_, nrow(g), ncol(g), dimnames=dimnames(g))
num[g=="AA"] <- 0L; num[g=="AB"] <- 1L; num[g=="BB"] <- 2L
write.table(num, "geno.tsv", sep="\t", quote=FALSE, na="NA", col.names=NA)
m <- cattle$map
write.table(data.frame(snp=rownames(m), chr=m$chr, pos_mb=m$pos), "map.tsv",
            sep="\t", quote=FALSE, row.names=FALSE)
ped <- cattle$pedigree
gped <- ped[ped$ID %in% rownames(g), ]
write.table(gped, "ped.tsv", sep="\t", quote=FALSE, row.names=FALSE)
RS
echo OK
