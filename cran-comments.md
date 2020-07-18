## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking installed package size ... NOTE
    installed size is  5.1Mb
    sub-directories of 1Mb or more:
      data   1.0Mb
      libs   3.4Mb

  The data in the sub-directory data/ has not changed since last version, while the size of the .so in the sub-directory libs/ is due to the compiled code linking against GSL.