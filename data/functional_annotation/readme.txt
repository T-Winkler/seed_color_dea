## Functional annotation

Functional annotation results from egg-NOG mapper, online tool usage 13.06.2022. Used the manually corrected polished genome annotation as input with the following settings:
## Mon Jun 13 13:33:07 2022
## emapper-2.1.7
## /data/shared/home/emapper/miniconda3/envs/eggnog-mapper-2.1/bin/emapper.py --cpu 20 --mp_start_method forkserver --data_dir /dev/shm/ -o out --output_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8wexw920 --temp_dir /emapper_web_jobs/emapper_jobs/user_data/MM_8wexw920 --override -m diamond --dmnd_ignore_warnings -i /emapper_web_jobs/emapper_jobs/user_data/MM_8wexw920/queries.fasta --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope 33090 --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel

# to enable correct loading of the file into R, I subset the relevant columns using the following command:
awk -F'\t' 'BEGIN {OFS = FS} {print $1,$10}' data/functional_annotation/MM_8wexw920.emapper.annotations.tsv > data/functional_annotation/egg_nog_mapper_go_terms.txt
