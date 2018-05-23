#!/usr/bin/env bash
cat /srv/data/genome/mouse/ensembl-92/Mus_musculus.GRCm38.92.gtf | grep 'protein_coding' > /home/xinhe/Projects/empirical_utrs/Mus_musculus.GRCm38.92.gtf
head -1000 /home/xinhe/Projects/empirical_utrs/Mus_musculus.GRCm38.92.gtf > /home/xinhe/Projects/empirical_utrs/Mus_musculus.GRCm38.92.test.gtf



python ./bin/get_empirical_utrs /srv/data/genome/mouse/ensembl-92/Mus_musculus.GRCm38.92.gtf \
/srv/data/eosterweil/sang_cage/hippocampus_adult.CNhs10478.13-16E8.mm10.nobarcode.bam \
--log-level debug &> log.txt



nohup python ./bin/get_empirical_utrs /home/xinhe/Projects/empirical_utrs/Mus_musculus.GRCm38.92.gtf /srv/data/eosterweil/sang_cage/hippocampus_adult.CNhs10478.13-16E8.mm10.nobarcode.bam /home/xinhe/Projects/empirical_utrs/output.txt --log-level debug >nohup.out.no.overlap &