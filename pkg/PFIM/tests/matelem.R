require(PFIM)

con <- file(system.file(package="PFIM", "Examples", "Example_5", "matelem.tmp"), open="r")
prot <- popProt$new(prots=readMatElem(con))
close(con)

prot$freq()                    # initialized to single protocol with maximum determinant
which(as.logical(prot$freq())) # 1-based index of protocol with maximum determinant
prot$props()                   # another way of getting the same

op <- options(digits=3)        # more concise printing of numeric values
prot$Fim()                     # Fisher Information matrix for these (trivial) proportions
prot$LMat()                    # Cholesky factor of the FIM
options(op)                    # restore the previous value

all.equal(tcrossprod(prot$LMat()), prot$Fim()) # check that LL' = F
prot$ldet()                    # log-determinant of the FIM
2 * log(prod(diag(prot$LMat()))) # another way of calculating the log-determinant of FIM

prot$ajout()                   # try to add another individual protocol but no other improves it

unzip(system.file(package="PFIM", "Examples", "ECo", "matelem1.zip"), "matelem1.tmp")
con <- file("matelem1.tmp", open="r")
prot <- popProt$new(prots=readMatElem(con))
close(con)
unlink("matelem1.tmp")

which(as.logical(prot$freq())) # 1-based index of protocol with maximum determinant
prot$ldet()                    # logarithm of that determinant

while (prot$ajout() && length(prot$props()) < 100L) {} # add up to 99 protocols

op <- options(digits=3)
rev(sort(prot$props()))        # combination (decreasing coefficients)
prot$Fim()                     # Fisher Information matrix
options(op)

prot$ldet()                    # log-determinant for FIM with 100 non-zeros

q('no')
