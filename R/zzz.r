# message so users don't get worried by the KernSmooth message when
# loading CorKern

cat("\nLoading GenKern version 0.98\n")
cat("Copyright Lucy and Aykroyd 2000\n")
cat("requires KernSmooth\n\n")

# required packages -  Wand and Jones' KernSmooth package mainly for 
# the dpik() function for the default h's
require(KernSmooth)

cat("\nPackage GenKern installed\n\n")

# load up the c module   
.First.lib <- function(lib, pkg) library.dynam("GenKern", pkg, lib)
# library.dynam("libCorKern")


