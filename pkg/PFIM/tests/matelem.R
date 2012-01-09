require(PFIM)

con <- file(system.file(package="PFIM", "Examples", "Example_5", "matelem.tmp"), open="r")
prot <- popProt$new(prots=readMatElem(con))
close(con)

prot$freq()                    # initialized to single protocol with maximum determinant
which(as.logical(prot$freq())) # 1-based index of protocol with maximum determinant
prot$ldet()                    # logarithm of that determinant
prot$Fim()                     # Fisher Information matrix for this (trivial) combination

prot$ajout()                   # try to add another but no other improves it

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
options(op)

prot$Fim()                     # Fisher Information matrix
prot$ldet()                    # log-determinant for FIM with 100 non-zeros

q('no')
