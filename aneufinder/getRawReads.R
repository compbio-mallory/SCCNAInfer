library(GenomicRanges)

# Function to process a single file
process_file <- function(file_path) {
  # Load the RData file
  load(file_path)
  print("load file")
  # Assuming the loaded object is named 'binned.data'
  gr <- bins[[1]]
  # Extract the required columns
  seqnames_list <- as.list(seqnames(gr))
  CHROM <- sapply(seqnames_list, as.character)
  START <- (start(gr))
  END <- (end(gr))
  counts <- mcols(gr)$counts
  # Get the sample name from the file name
  file_name<- gsub("^.*?([^/]*)_.*$", "\\1", file_path)
  sample_name <- strsplit(file_name, "\\.")[[1]][1]
  # Combine the extracted columns into a data frame
  result <- data.frame(CHROM, START, END)
  result[sample_name] <- counts
  print(paste("Got result for", sample_name))
  return(result)
}

# Specify the directory containing the RData files
# Get the directory path from the command-line arguments
args <- commandArgs(trailingOnly = TRUE)
directory <- args[1]

# Initialize an empty data frame for the final result
final_result <- data.frame()

# Iterate through each file in the directory
for (file_name in list.files(directory, pattern = "\\.RData$", full.names = TRUE)) {
  # Process the current file
  print(file_name)
  current_result <- process_file(file_name)
  
  # Merge the current result with the final result
  if (nrow(final_result) == 0) {
    final_result <- current_result
  } else {
    final_result <- merge(final_result, current_result, by = c("CHROM", "START", "END"), all = TRUE)
  }
}
sorted_df <- final_result[order(final_result$CHROM, final_result$START), ]
write.table(sorted_df, args[2], sep = "\t", row.names = FALSE, quote = FALSE)

