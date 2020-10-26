library(ggplot2)

FPKMVsCountsPlot <- function(fpkm.merged.filename, counts.merged.filename, description, eq.place.x, eq.place.y, regression = "linear", sample.limit = 1000) {
  fpkm.table <- read.table(file = fpkm.merged.filename,
                           row.names = 1,header = TRUE,sep = "\t",
                           colClasses = "character", check.names = FALSE,comment.char = "",
                           fill = TRUE, blank.lines.skip = FALSE)
  counts.table <- read.table(file = counts.merged.filename,
                             row.names = 1,header = TRUE,sep = "\t",
                             colClasses = "character", check.names = FALSE,comment.char = "",
                             fill = TRUE, blank.lines.skip = FALSE)
  counts <- c()
  fpkm <- c()
  for (r in row.names(fpkm.table)[1:sample.limit]) {
    counts <- c(counts, as.numeric(counts.table[r, ]))
    fpkm <- c(fpkm, as.numeric(fpkm.table[r, ]))
  }
  if (regression == "linear") {
    regression.line <- lm(fpkm ~ counts)
    coefficients <- coef(regression.line)
    plot <- qplot(x=counts, y = fpkm) +
            geom_abline(intercept = coefficients[1],
                        slope = coefficients[2]) +
            annotate("text", x = eq.place.x, y = eq.place.y,
                     label = lm_eqn(counts, fpkm), parse = T)
  } else if (regression == "logarithmic") {
    fit <- nls(fpkm ~ a + b * log(counts + 1), start = list(a = 1, b = 1))
    p <- qplot(x=counts, y=fpkm) +
            annotate("text", x = eq.place.x, y = eq.place.y,
                     label = log_eqn(counts, fpkm), parse = T) +
            stat_smooth(formula =
                        #FIX HERE
    matlines(x=seq(from=1,to=max(counts[1:sample.limit]),length.out=1000),
             y=predict(fit,newdata=list(counts=seq(from=1,to=max(counts[1:sample.limit]),length.out=1000))))
  }
  ggsave(plot, file = paste0("fpkm_vs_counts.", description, ".jpg"))
  print(cor(counts, fpkm, use = "complete.obs"))
}

lm_eqn <- function(x, y){
  m <- lm(y ~ x);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq))
}

log_eqn <- function(x, y){
  m <- nls(y ~ a + b * log(x + 1), start = list(a = 1, b = 1))
  RSS <- sum(residuals(m)^2)
  TSS <- sum((y - mean(y))^2)
  r.squared <- 1 - (RSS/TSS)
  eq <- substitute(italic(y) == a + b %.% italic(ln(x))*","~~italic(r)^2~"="~r2,
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),
                        r2 = format(r.squared, digits = 3)))
  as.character(as.expression(eq))
}
