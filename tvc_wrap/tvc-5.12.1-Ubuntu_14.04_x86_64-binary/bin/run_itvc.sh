python ./variant_caller_pipeline.py \
	-b /home/aloe/src/genetic_playground/data/CCP.designed.bed \
	-i /home/aloe/src/genetic_playground/data/mil-2020-CCP-1375-pcblt.bam \
	-r /home/aloe/src/genetic_playground/data/hg19_25.fa \
	-o /home/aloe/src/genetic_playground/tvc_out \
	-p /home/aloe/src/genetic_playground/data/CCP.ITVC.json \
	-N 15