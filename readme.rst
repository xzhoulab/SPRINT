
##SPRINT

**SPRINT** is an efficient method to identify genes with spatial expression pattern. 
The intended applications are spatially resolved RNA-sequencing from e.g.
Spatial Transcriptomics, or *in situ* gene expression measurements from
e.g. SeqFISH, Merfish.

## System requirements
SPRINT has been tested on R 3.3.1 and is platform independent (tested on Linux, OS X and Windows)
Installation can then be done via the devtools package:

```R
library('devtools')
devtools::install_github('xzlab/SPRINT')
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

    .. raw:: html

    <div>
    <style scoped>
        .dataframe tbody tr th:only-of-type {
            vertical-align: middle;
        }
    
        .dataframe tbody tr th {
            vertical-align: top;
        }
    
        .dataframe thead th {
            text-align: right;
        }
    </style>
    <table border="1" class="dataframe">
      <thead>
        <tr style="text-align: right;">
          <th></th>
          <th>GAPDH</th>
          <th>USP4</th>
          <th>MAPKAPK2</th>
          <th>CPEB1</th>
          <th>LANCL2</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <th>17.907x4.967</th>
          <td>1</td>
          <td>1</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>18.965x5.003</th>
          <td>7</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>18.954x5.995</th>
          <td>5</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>17.846x5.993</th>
          <td>1</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
          <td>0</td>
        </tr>
        <tr>
          <th>20.016x6.019</th>
          <td>2</td>
          <td>0</td>
          <td>1</td>
          <td>0</td>
          <td>0</td>
        </tr>
      </tbody>
    </table>
    </div>

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










