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
  select(-singlem_fraction, -average_genome_size)

new_smf <- read_delim("singlem_new.tsv", col_names = c("EHI_number", "singlem_fraction", "average_genome_size"))

#Now merge the tables with the correct values, then update AirTable
merged <- airtable_data %>%
  inner_join(., new_smf, by = join_by(EHI_plaintext == EHI_number)) %>%
  #function for updating records in AirTable
  update_records(airtable = table, 
                 airtable_id_col = airtable_record_id,
                 columns = c(singlem_fraction, average_genome_size))

#Create plot comparing old to new SMF values
figure <- read_delim("singlem_input.csv") %>%
  inner_join(new_smf, by = join_by(alias_fastq == EHI_number)) %>%
  ggplot(aes(x = singlem_fraction_before, y = singlem_fraction)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  ylab("singlem_fraction_after")

ggsave("comparison.png", figure, width = 8, height = 8)
