mytextplot <- function(x, y, words, cex = 1, new = TRUE, show.lines = TRUE,
                        ...){
    if (new)
        plot(x, y, type = "n", ...)
    lay <- wordlayout(x, y, words, cex, ...)
    if (show.lines) {
        for (i in 1:length(x)) {
            xl <- lay[i, 1]
            yl <- lay[i, 2]
            w <- lay[i, 3]
            h <- lay[i, 4]
            if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] >
                yl + h) {
#                points(x[i], y[i], pch = 16, col = "red", cex = 0.5)
                nx <- xl + 0.5 * w
                ny <- yl + 0.5 * h
                lines(c(x[i], nx), c(y[i], ny), col = "grey")
            }
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3], lay[, 2] + 0.5 * lay[, 4],
         words, cex = cex, ...)
}
