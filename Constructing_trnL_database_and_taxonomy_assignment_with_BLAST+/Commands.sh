# Download sequences
esearch -db nucleotide -query "biomol_genomic[PROP] AND is_nuccore[filter] AND chloroplast[filter] AND (asteraceae OR brassicaceae OR chenopodiaceae OR fabaceae OR plantaginaceae OR poaceae OR cupressaceae OR zygophyllaceae)" | efetch -db nucleotide -format fasta -mode xml > Sequences.Raw.2020.06.03.xml

# Modify the XML so it can be parsed in Python
head -3 Sequences.Raw.2020.06.03.xml > Sequences.2020.06.03.xml
grep -v '<?xml\|<!DOCTYPE\|TSeqSet' Sequences.Raw.2020.06.03.xml >> Sequences.2020.06.03.xml
tail -1 Sequences.Raw.2020.06.03.xml >> Sequences.2020.06.03.xml

# Download taxonomies; there are 32,661 unique taxonomy IDs, and efetch downloads 10,000 at a time by default, so download them in 4 pieces
efetch -db taxonomy -id "$(grep 'TSeq_taxid' Sequences.2020.06.03.xml | cut -d '>' -f 2 | cut -d '<' -f 1 | sort -u | head -10000 | tr '\n' ';')" -format xml > Taxonomies.01.2020.06.03.xml
efetch -db taxonomy -id "$(grep 'TSeq_taxid' Sequences.2020.06.03.xml | cut -d '>' -f 2 | cut -d '<' -f 1 | sort -u | head -20000 | tail -10000 | tr '\n' ';')" -format xml > Taxonomies.02.2020.06.03.xml
efetch -db taxonomy -id "$(grep 'TSeq_taxid' Sequences.2020.06.03.xml | cut -d '>' -f 2 | cut -d '<' -f 1 | sort -u | head -30000 | tail -10000 | tr '\n' ';')" -format xml > Taxonomies.03.2020.06.03.xml
efetch -db taxonomy -id "$(grep 'TSeq_taxid' Sequences.2020.06.03.xml | cut -d '>' -f 2 | cut -d '<' -f 1 | sort -u | head -32661 | tail -2661 | tr '\n' ';')" -format xml > Taxonomies.04.2020.06.03.xml

# Generate a combined FASTA file
Scripts/xmls_to_fasta.py Sequences.2020.06.03.xml genus environmental_samples Taxonomies.0*.xml > Combined.2020.06.03.fasta

# Trim the combined FASTA file with Cutadapt 2.10
cutadapt -g "GGGCAATCCTGAGCCAA;max_error_rate=0.12" -a "GATAGGTGCAGAGACTCAATGG;max_error_rate=0.1" --no-indels -n 2 --action=lowercase --discard-untrimmed Combined.2020.06.03.fasta > Cutadapt.2020.06.03.fasta

# Download Cactaceae sequences
esearch -db nucleotide -query "biomol_genomic[PROP] AND is_nuccore[filter] AND chloroplast[filter] AND cactaceae" | efetch -db nucleotide -format fasta -mode xml > Sequences.Cactaceae.Raw.2021.09.29.xml

# Download Cactaceae taxonomies; there are 1,760 unique taxonomy IDs, and efetch downloads 10,000 at a time by default, so download them all at once
efetch -db taxonomy -id "$(grep 'TSeq_taxid' Sequences.Cactaceae.2021.09.29.xml | cut -d '>' -f 2 | cut -d '<' -f 1 | sort -u | tr '\n' ';')" -format xml > Taxonomies.Cactaceae.2021.09.29.xml

# Modify the Cactaceae XML so it can be parsed in Python
head -3 Sequences.Cactaceae.Raw.2021.09.29.xml > Sequences.Cactaceae.2021.09.29.xml
grep -v '<?xml\|<!DOCTYPE\|TSeqSet' Sequences.Cactaceae.Raw.2021.09.29.xml >> Sequences.Cactaceae.2021.09.29.xml
tail -1 Sequences.Cactaceae.Raw.2021.09.29.xml >> Sequences.Cactaceae.2021.09.29.xml

# Generate a combined FASTA file for the Cactaceae sequences
Scripts/xmls_to_fasta.py Sequences.Cactaceae.2021.09.29.xml genus environmental_samples Taxonomies.Cactaceae.2021.09.29.xml > Combined.Cactaceae.2021.09.29.fasta

# Trim the combined Cactaceae FASTA file with Cutadapt 2.10
cutadapt -g "GGGCAATCCTGAGCCAA;max_error_rate=0.12" -a "GATAGGTGCAGAGACTCAATGG;max_error_rate=0.1" --no-indels -n 2 --action=lowercase --discard-untrimmed Combined.Cactaceae.2021.09.29.fasta > Cutadapt.Cactaceae.2021.09.29.fasta

# Build a database FASTA using the two Cutadapt-trimmed FASTAs
Scripts/build_db.py <(cat Cutadapt.2020.06.03.fasta Cutadapt.Cactaceae.2021.09.29.fasta) GGGCAATCCTGAGCCAA GATAGGTGCAGAGACTCAATGG 175 8 > Database.2021.09.29.fasta

# Build a BLAST database so the OTU FASTA sequences can be aligned
makeblastdb -in Database.2021.09.29.fasta -dbtype nucl

# Convert the OTU tables to FASTAs
Scripts/csv_to_fasta.py creosote-juniper-OTU-table.csv > creosote-juniper-OTU.fasta
Scripts/csv_to_fasta.py <(sed 's/"//g' cactus-OTU-table.csv) > cactus-OTU.fasta

# Since some of the sequences are under 30 nucleotides, align using the "blastn-short" task with a maximum E-value of 1 and a percent identity of 90
blastn -task "blastn-short" -db Database.2021.09.29.fasta -query creosote-juniper-OTU.fasta -evalue 1 -word_size 4 -perc_identity 90 -outfmt 6 > creosote-juniper-OTU.outfmt6
blastn -task "blastn-short" -db Database.2021.09.29.fasta -query cactus-OTU.fasta -evalue 1 -word_size 4 -perc_identity 90 -outfmt 6 > cactus-OTU.outfmt6

# Assign taxonomy based on the BLAST results, requiring the hit and OTU sequence lengths to be within 10% of each other, and taxonomy to be assigned at a given rank when 50% or more of the BLAST hits agree
Scripts/assign_taxonomy.py creosote-juniper-OTU.outfmt6 Database.2021.09.29.fasta creosote-juniper-OTU.fasta 0.9 0.5 > creosote-juniper-OTU.Assigned.csv
Scripts/assign_taxonomy.py cactus-OTU.outfmt6 Database.2021.09.29.fasta cactus-OTU.fasta 0.9 0.5 > cactus-OTU.Assigned.csv
