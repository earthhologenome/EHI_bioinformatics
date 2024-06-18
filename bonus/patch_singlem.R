#Script for patching new SingleM microbial fraction values to AirTable
#Raphael Eisenhofer 2024

library(readr)
library(magrittr)
library(dplyr)
library(ggplot2)
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
  select(-singlem_fraction, -average_genome_size)

new_smf <- read.delim("singlem_new.tsv", sep = "\t", header = FALSE, 
                      col.names = c("EHI_number", "singlem_fraction", "average_genome_size"))


tryCatch({
  result <- airtable_data %>%
    inner_join(., new_smf, by = join_by(EHI_plaintext == EHI_number)) %>%
    #function for updating records in AirTable
    update_records(airtable = view, 
                   airtable_id_col = airtable_record_id,
                   columns = c(singlem_fraction, average_genome_size),
                   safely = FALSE)
  print(result)
}, error = function(e) {
  # Print the error message
  print(e)
  
  # Additional debugging
  if (inherits(e, "HTTPError")) {
    response <- e$response
    print(content(response, "text", encoding = "UTF-8"))
  }
})


#Create plot comparing old to new SMF values
figure <- read.delim("singlem_input.csv", sep = ",") %>%
  inner_join(new_smf, by = join_by(alias_fastq == EHI_number)) %>%
  ggplot(aes(x = singlem_fraction_before, y = singlem_fraction)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  ylab("singlem_fraction_after")

ggsave("comparison.png", figure, width = 8, height = 8)
