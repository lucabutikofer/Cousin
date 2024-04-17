
# Function to clean species names in a list of species names to conform them to the foollowing format:
# 
# 1. Genus speices
# 2. Genus species subspecies
# 
# No other words, acronyms and formats will be left.
# 
# require(dplyr)


clean.sp.names <- function(x){
  
  # x: Character vector of species names
  
  
  # remove all brackets and their content (e.g. "(L.)")
  x <- gsub("\\s*\\([^\\)]+\\)", "", x)
  
  # remove all words after the one following "subsp." (typically the authority)
  x <-
    strsplit(x, " ") %>%
      sapply(function(X){
        if (any(grepl("subsp.", X))){ # if string contains "subsp."
          if (length(X) > (which(grepl("subsp.", X)) + 1)){ # if there are words after the ssp. designation
            X <- X[1:(which(grepl("subsp.", X)) + 1)]
          } else X
        } else X
        X <- paste(X, collapse = " ")
      }) %>% gsub("  ", " ", .)

  # remove "subsp."
  x <- gsub("subsp. ", "", x)
  
  # remobe all characters between a space and a "." (e.g. "Linn." or "L.", typically authority abbreviations)
  x <- strsplit(x, " ") %>% sapply(function(X) X[!grepl("\\.", X)] %>% paste(collapse = " "))
  
  # remove words with capitalised first letter that are not the first word (e.g. "Linnaeus" but not "Pinus")
  x <-
    strsplit(x, " ") %>%
    sapply(function(X){
      retain <- which(!grepl("[A-Z]", X))
      retain <- c(1, retain)
      X <- X[retain]
      return(paste(X, collapse = " "))
    })
  
  # remove "&"
  x <- gsub("&", "", x)
  
  # remove double and trailing spaces
  x <- gsub("  ", " ", x)
  x <- gsub(" $", "", x)
  
  # replace "×" with "x" (yes, they are different!)
  x <- gsub("×", "x", x)
  
  # attach " x " to specific epithet (e.g. Salix ×fragilis instead of Salix x fragilis)
  x <- gsub(" x ", " x", x)
  
  return(x)
  
}
