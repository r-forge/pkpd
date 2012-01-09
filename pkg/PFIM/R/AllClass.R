### Class definitions for the PFIM package

popProt <-
    setRefClass("popProt",
                fields  = list(prots = "list", Ptr = "externalptr"),
                methods =
                list(
                     ajout = function() {
                         'One additive step in the Fedorov-Wynn algorithm'
                         .Call(PopProt_ajout, ptr())
                     },
                     Fim = function() {
                         'returns the Fisher information matrix'
                         .Call(PopProt_fim, ptr())
                     },
                     freq = function() {
                         'returns the vector of frequencies'
                         .Call(PopProt_freq, ptr())
                     },
                     iTime = function(i, base1 = TRUE) {
                         'returns the vector of times for the ith protocol'
                         .Call(PopProt_iTime, ptr(), as.integer(i) - as.integer(base1))
                     },
                     iFim = function(i, base1 = TRUE) {
                         'returns the Fisher information matrix for the ith protocol'
                         .Call(PopProt_iFim, ptr(), as.integer(i) - as.integer(base1))
                     },
                     LMat = function() {
                         'returns the Cholesky factor of the Fisher information matrix'
                         .Call(PopProt_Lmat, ptr())
                     },
                     ldet = function() {
                         'returns the logarithm of the determinant of the Fim'
                         .Call(PopProt_ldet, ptr())
                     },
                     normalize = function() {
                         'normalize the proportions and drop those below the threshold'
                         .Call(PopProt_normalize, ptr())
                     },
                     nprot = function() {
                         'returns the number of protocols available'
                         .Call(PopProt_nprot, ptr())
                     },
                     props = function() {
                         'return the nonzero proportions'
                         prop <- .Call(PopProt_freq, ptr())
                         names(prop) <- seq_along(prop)
                         prop[which(prop > 0)]
                     },
                     ptr          = function() {
                         'returns the external pointer, regenerating if necessary'
                         if (length(prots)) {
                             if (.Call(isNullExtPtr, Ptr))
                                 Ptr <<- .Call(PopProt_Create, prots);
                         }
                         Ptr
                     },
                     setThreshold = function(thresh) {
                         'set the threshold for declaring a proportion to be zero'
                         .Call(PopProt_setThreshold, ptr(), thresh)
                     },
                     threshold = function() {
                         'the threshold for declaring a proportion to be zero'
                         .Call(PopProt_threshold, ptr())
                     })
                )

popProt$lock("prots")
