## Following comands are used to analyze the chimeras after running the vserach pipline for OTU clustering.
## Kabin Xie
## 2020.5.4

## all.preculstered.fasta are generated in VSEARCH OTU clustering pipeline.
#detect chimera using uchime_denovo
vsearch  --uchime_denovo all.preclustered.fasta --sizein --sizeout --fasta_width 0 --nonchimeras all.denovo.nonchimeras.2.fasta --chimeras  all.denovo.chimeras.fasta  --uchimeout all.denovo.chimeras.uchiout

# Extract all chimeric sequences, dereplicated
../map.pl all.derep.fasta all.preclustered.uc all.denovo.chimeras.fasta > all.denovo.chimeras.derep.fasta

#Extract all chimeric sequences in each samples
../map.pl all.fasta all.derep.uc all.denovo.chimeras.derep.fasta > all.chimeras.fasta

# Generate the chimera table (same as the OTU table)
vsearch --cluster_size all.chimeras.fasta \
    --id 1 \
    --strand plus \
    --minsize 8 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.clustered.uc \
    --relabel Chi_ \
    --centroids all.chimeras100.fasta \
    --otutabout all.chimeras100.txt

##Done!
