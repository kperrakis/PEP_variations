reset <-
  function() {
    par(
      mfrow = c(1, 1),
      oma = rep(0, 4),
      mar = rep(0, 4),
      new = TRUE
    )
    plot(
      0:1,
      0:1,
      type = "n",
      xlab = "",
      ylab = "",
      axes = FALSE
    )
  }
