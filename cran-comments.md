## Test environments
* local MacOS 10.13 install, R 3.4.3
* MacOS (on travis-ci), R 3.3.3 and 3.4.4
* Ubuntu 14.04 (on travis-ci), R 3.3.3 and 3.4.4

## R CMD check results
 - MacOS (local and travis-ci): There were no ERRORs, WARNINGs, or NOTEs. 

 - Ubuntu (travis-ci): There were no ERRORs or WARNINGs.
There was 1 NOTE:
| checking installed package size ... NOTE
|   installed size is  5.7Mb
|   sub-directories of 1Mb or more:
|     libs   5.6Mb

Try with devtools::check()
