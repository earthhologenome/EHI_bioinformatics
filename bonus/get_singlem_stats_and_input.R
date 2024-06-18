#Script for pulling current SingleM values from AirTable and for creating the 
#input file for update_singlem_airtable.snakefile
#Raphael Eisenhofer 2024
library(magrittr)
library(dplyr)
library(rairtable)

#setup api key
apikey <- read_file("/projects/ehi/data/.airtable_api_key.txt")
set_airtable_api_key(apikey)

#Select the right view id
##Note that I've filtered out samples that give NaN or #ERROR in formulas, as 
##these prevent rAirtable from pulling the data. See 'singlem_update' view
view <- airtable('tblJfLRU2FIVz37Y1', 
                 'appQpr6MxnaiVHsHy', 
                 view = 'viwYBdUy5zN7ZFAMG')

airtable_data <- read_airtable(view, id_to_col = TRUE, max_rows = 50000) %>%
  select(alias_fastq, singlem_fraction, average_genome_size, URL_meta1, URL_meta2) %>%
  rename(singlem_fraction_before = singlem_fraction, average_genome_size_before = average_genome_size) %>%
  mutate(alias_fastq = str_replace(alias_fastq, "M.", ""))

#save output
write.csv(x = airtable_data, file = "singlem_input.csv", row.names = F)
