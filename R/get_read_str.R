get_read_str <- function(protocol) {
  original_protocol <- protocol
  error_msg <- paste("No read structure for", original_protocol, "found")

  # strip non-alphanumeric characters and convert to lower case
  protocol <- gsub("[^[:alnum:]_]", "", protocol)
  protocol <- tolower(protocol)

  switch(
    protocol,
    "celseq" = list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8, us = 0, ul = 6),
    "celseq2" = list(bs1 = -1, bl1 = 0, bs2 = 6, bl2 = 8, us = 0, ul = 6),
    "dropseq" = list(bs1 = -1, bl1 = 0, bs2 = 0, bl2 = 12, us = 12, ul = 8),
    default = stop(error_msg)
  )
}

