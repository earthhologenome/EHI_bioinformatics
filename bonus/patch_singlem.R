#Script for patching new SingleM microbial fraction values to AirTable
#Raphael Eisenhofer 2024

library(tidyverse)
library(rairtable)

#Select the right view id
##Note that I've filtered out samples that give NaN or #ERROR in formulas, as 
##these prevent rAirtable from pulling the data. See 'singlem_update' view
view <- airtable('tblJfLRU2FIVz37Y1', 
                 'appQpr6MxnaiVHsHy', 
                 view = 'viwYBdUy5zN7ZFAMG')

airtable_data <- read_airtable(view, id_to_col = TRUE, max_rows = 50000) %>%
  select(alias_fastq, singlem_fraction, URL_meta1, URL_meta2) %>%
  rename(singlem_fraction_before = singlem_fraction) %>%
  mutate(alias_fastq = str_replace(alias_fastq, "M.", ""))

#save output
write.csv(x = airtable_data, file = "singlem_input.csv", row.names = F)