
## SPRINT

**SPRINT** is an efficient method to identify genes with spatial expression pattern. 
The intended applications are spatially resolved RNA-sequencing from e.g.
Spatial Transcriptomics, or *in situ* gene expression measurements from
e.g. SeqFISH, Merfish.

## System requirements
SPRINT has been tested on R 3.3.1 and is platform independent (tested on Linux, OS X and Windows)
Installation can then be done via the devtools package:

```R
library('devtools')
devtools::install_github('xzhoulab/SPRINT')
```
Alternatively, installation can then be done from a local binary package
tarball from the shell:
```bash
R CMD INSTALL SPRINT_1.0.0.tar.gz
```



## Sample Code: Analysis of Breast Cancer Data
```R
    library('SPRINT')
    load("~/data/Layer2_BC_Count.rds")
     
    ## rawcount matrix of genes by cells/spots
    rawcount[1:5,1:5]
    
    ## extract the coordinates from the rawdata
    info <- cbind.data.frame(x=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",1)),
                            y=as.numeric(sapply(strsplit(colnames(rawcount),split="x"),"[",2)),
                            total_counts=apply(rawcount,2,sum))
    rownames(info) <- colnames(rawcount)

    ## filter genes and cells/spots and create the SPRINT object for following analysis
    sprint <- CreateSPRINTObject(counts=rawcount, location=info[,1:2])

    ## total counts for each cell/spot
    sprint@lib_size <- apply(sprint@counts, 2, sum)

    ## Take the first ten genes as an example
    sprint@counts   <- sprint@counts[1:10,]

    ## Estimating Parameter Under Null
    sprint <- sprint.vc(sprint,covariates=NULL, lib_size=sprint@lib_size, num_core=1,verbose=F)

    ## Calculating pval
    sprint <- sprint.test(sprint, check_positive = T, verbose=F)
```










