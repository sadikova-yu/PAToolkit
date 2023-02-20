use strict;
use warnings;

print `bedtools intersect -b $ARGV[0] -a /home/aod/RnD/somaticRecall/CosmicCodingMuts.chr.vcf | grep -v "ENST" | grep -v "ENSG" | awk '{if (length(\$4)<30 && length(\$5) < 30) print}' | sort -k1,1 -k2,2n | uniq `;
