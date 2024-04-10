# helper function for plotting troll world object as trials unfold
# requires that ffmpeg be on the system path
plot_troll_world <- function(tw, width = 8, height = 5, ncores = 8, frame_rate = 5, scale_y = "fixed", out_mp4 = "trolls.mp4") {
  require(magick)
  require(ggplot2)
  require(glue)
  require(parallel)
  require(foreach)
  require(checkmate)

  checkmate::assert_class(tw, "troll_world")

  # combine aspects of the contingency for plotting
  vv <- tw$get_values_matrix()
  vdf <- reshape2::melt(vv, varnames = c("trial", "point"))
  ee <- tw$erasure_segments
  pos_df <- data.frame(pos_rad = tw$get_pvec(), point = seq_along(tw$get_pvec()))

  drift_df <- data.frame(trial=seq_along(tw$get_drift_vec()), drift = tw$get_drift_vec())
  drift_df$drift_rad <- drift_df$drift*diff(tw$get_pvec()[1:2]) # convert to radians based on d_theta between points

  vdf <- vdf %>%
    left_join(ee) %>%
    left_join(pos_df) %>%
    left_join(drift_df)  
  
  tvec <- sort(unique(vdf$trial))
  images <- rep(NA_character_, length(tvec))

  if (ncores > 1L) {
    cl <- makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    on.exit(try(stopCluster(cl)))
  } else {
    registerDoSEQ()
  }

  if (scale_y == "fixed") {
    lims <- c(min(vdf$value), max(vdf$value))
  } else {
    lims <- NULL
  }
  minval <- min(vdf$value)

  images <- foreach(ii = tvec, .combine = c, .packages = c("dplyr", "ggplot2", "glue")) %dopar% {
    dd <- vdf %>% filter(trial == {{ ii }})
    this_cond <- dd$trial_type[1]
    ss <- ifelse(this_cond != "no erasure", glue("min: {round(dd$segment_min_shift[1], 1)}, max: {round(dd$segment_max_shift[1], 1)}"), "")
    ss <- ifelse(this_cond == "erasure", glue("{ss}, pct_change: {round((dd$v_new[1] - dd$v_old[1])/abs(dd$v_old[1]) * 100, 2)}"), ss) # add percent change for erasures

    gg <- ggplot(dd, aes(x = pos_rad, y = value)) +
      geom_line() +
      ggtitle(glue("Trial: {ii}"), subtitle = glue("Condition: {this_cond} {ss}, drift: {dd$drift[1]}, drift_rad: {round(dd$drift_rad[1], 2)}")) +
      ylim(lims)

    if (this_cond != "no erasure") {
      cc <- ifelse(this_cond == "attention", "blue", "orange")
      ymax <- ifelse(is.null(lims), max(dd$value), lims[2]) # rectangle encompasses y range
      rect_df <- dd[1, ]
      if (rect_df$segment_max_shift - rect_df$segment_min_shift < 0) { # handle wraparound at 0
        rect_df <- rbind(rect_df, rect_df) # make 2 rows to draw segments near 0 and 2pi
        rect_df$segment_min_shift[1L] <- 0 # first row is 0 -> low
        rect_df$segment_max_shift[2L] <- 2 * pi # second row is high -> 2pi
      }
      gg <- gg +
        geom_rect(
          data = rect_df,
          mapping = aes(xmin = segment_min_shift, xmax = segment_max_shift, ymin = minval, ymax = ymax),
          fill = cc, color = "transparent", alpha = 0.2
        )
    }

    tmpfile <- tempfile(pattern = paste0("img", sprintf("%05d", ii), "_"), fileext = ".png")

    ggsave(
      filename = tmpfile,
      plot = gg,
      device = "png",
      width = width, height = height,
      units = "in"
    )

    return(tmpfile)
  }

  # write image list to file
  flist <- tempfile(fileext = ".txt")

  # need duration for each image if using a flat file: https://stackoverflow.com/questions/76988814/ffmpeg-combine-images-to-video-and-specify-duration-for-each-image
  writeLines(paste(glue("file {images}"), glue("duration {1/frame_rate}"), sep = "\n"), con = flist)

  # ffmpeg -f concat -i <(for f in "$PWD"/${subid}*.MP4; do echo "file '$f'"; done) -c copy "$concatFile"
  # system(glue("ffmpeg -y -f concat -i <(for f in {flist}; do echo \"file '$f'\"; done) -framerate 1  -r 5 -c:v libx264 -pix_fmt yuv420p ~/video.mp4"))
  system2("ffmpeg", glue("-y -f concat -safe 0 -i {flist} -c:v libx264 -pix_fmt yuv420p {out_mp4}"), stdout = NULL, stderr = NULL) #-threads 24

  # clean up temp files
  on.exit(try(unlink(c(images, flist))), add = TRUE)

  # this is rather slow and not amenable to parallelization due to external pointers used by image_read
  # frames <- foreach(i = length(images):1, .packages="magick") %dopar% {
  # frames <- c()
  # for (i in length(images):1) {
  #  x <- image_read(images[i])
  #  x <- image_scale(x, "300")
  #  frames <- c(x, frames)
  # }

  # animation <- image_animate(frames, fps = 4)
  # image_write(animation, "plot.gif")

  # Use image magick (much slower than ffmpeg)
  # system("convert -delay 80 *.png animated_count_down.gif")
}
