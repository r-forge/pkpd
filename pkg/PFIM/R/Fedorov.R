##' Read a \dQuote{matelem} file written by the PFIM script.
##'
##' @title Read a \dQuote{matelem} file
##' @param con file name or connection
##' @return a list of lists.  The inner lists are elementary protocols.
readMatElem <- function(con) {
    if (is.character(con)) {
        con <- file(con, "r")
        on.exit(close(con))
    }
    ans <- list()
    el <- 1L
    Fd <- 0L
    while (length(ll <- readLines(con, n = 1L))) {
        if (scan(textConnection(ll), what=integer(), n=1L, quiet=TRUE) != el)
            stop("Format of \"", file, "\" seems corrupt")
        npts <- scan(con, integer(), n=1L, quiet=TRUE)
        tpts <- scan(con, double(), n=npts, quiet=TRUE)
        if (!Fd) {
            Fd <- length(r1 <- scan(textConnection(readLines(con, n=1L)), what=double(), quiet=TRUE))
        } else {
            r1 <- scan(con, double(), n=Fd, quiet=TRUE)
        }
        ans[[el]] <- list(times = tpts,
                          Fim = matrix(c(r1, scan(con, double(), n=Fd*(Fd-1L), quiet=TRUE)),
                          nc=Fd, byrow=TRUE))
        el <- el + 1L
    }
    ans
}

