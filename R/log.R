
log.end_time <- function(log_file, start_time, end_time) {

  capture.output(cat(paste("Script successfully finished at:", end_time, sep = " ")), file = log_file, append = T)

  capture.output(cat("\n\nTotal time:", round(difftime(end_time, start_time, units = "min"), 2), "minutes", sep = " "),
                 file = log_file, append = T)

}



log.options_summary <- function(log_file, start_time, gms, mds, opt, cl = FALSE){

  if(cl == TRUE){

    capture.output(cat(paste(strsplit(commandArgs(trailingOnly = FALSE)[4],"=")[[1]][2],
                             "start time:", start_time, "\n\n", sep = " ")),
                   file = log_file, append = T)

  } else {

    capture.output(cat(paste("start time:", start_time, "\n\n", sep = " ")),
                   file = log_file, append = T)

  }

  capture.output(cat("Gene Matrix Files:\n\n"), file = log_file, append = T)

  capture.output(cat(unlist(gms), sep = "\n"), file = log_file, append = T)

  capture.output(cat("\nResults Files:\n\n"), file = log_file, append = T)

  capture.output(cat(unlist(mds), sep = "\n"), file = log_file, append = T)

  opt[sapply(opt, hanta:::is.empty)] <- ""

  opt <- lapply(opt, toString)

  t <- t(data.frame(opt))

  colnames(t) <- "options"

  capture.output(cat("\n"), file = log_file, append = T)

  capture.output(t, file = log_file, append = T)

  capture.output(cat("\n"), file = log_file, append = T)

}



log.out <- function(body, time = TRUE, n = 1, file, append, echo = TRUE){

  if(time == TRUE){

    body = paste(body, Sys.time(), sep = " ")

  }

  body = paste0(body, strrep("\n", n))

  capture.output(cat(paste0(body)), file = file, append = append)

  if(echo == TRUE){

    write(body, stdout())

  }

}
