cdcat 1.0.6:

## cdcatR 1.0.6 release (May 25, 2022)
* The {cdcatR} package was removed from CRAN due to an issue with the {CDM} package. The {CDM
package has been successfully re-submitted to CRAN so there should not be any problem now
* Following the indications received on May 24, we now use inherits() instead of if() conditions comparing class() to string. Thank you. Please find below an explanation for the other two notes we get in some checks

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results
* There were no ERRORs or WARNINGs.
* There were 2 NOTES:
* First note: 
```
Found the following files/directories: 'lastMiKTeXException'
```
Please note that this might be due to a bug/crash in MiKTeX, as is noted in R-hub issue #503, which can be found at r-hub/rhub#503

* Second note:
```
Maintainer: 'Miguel A. Sorrel <miguel.sorrel@uam.es>'

  New submission
  
  Package was archived on CRAN
  
  Uses the superseded package: 'doSNOW'
  
   Found the following (possibly) invalid DOIs:
     DOI: 10.1177/0146621614554650
       From: DESCRIPTION
       Status: Service Unavailable
     DOI: 10.1177/0146621618813113
       Status: Service Unavailable
       Message: 503
       From: DESCRIPTION
```
Please note these two DOIs are valid, both are from the journal Applied Psychological Measurement. As indicated the package was archived due to a problem with the CDM package but that package made it to CRAN again successfully on May 13th. Finally the `doSNOW` package is required to show a progress bar inside a foreach loop

## cdcatR 1.0.5 release (April 5, 2022)
* We included a new item selection methods and conducted some minor modifications to adapt the outputs
* I am aware that I just updated the package a few days ago but I am teaching a course in the next few days and would like to have this new method. I have notified Prof. Uwe Ligges (thanks for your response!). I will be more careful in the future with more spaced out updates

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results
* There were no ERRORs or WARNINGs.
* There were 2 NOTES:
* First note: 
```
Found the following files/directories: 'lastMiKTeXException'
```
Please note that this might be due to a bug/crash in MiKTeX, as is noted in R-hub issue #503, which can be found at r-hub/rhub#503

* Second note:
```
Maintainer: 'Miguel A. Sorrel <miguel.sorrel@uam.es>'

  Days since last update: 5
  
  Uses the superseded package: 'doSNOW'
  
   Found the following (possibly) invalid DOIs:
     DOI: 10.1177/0146621614554650
       From: DESCRIPTION
       Status: Service Unavailable
     DOI: 10.1177/0146621618813113
       Status: Service Unavailable
       Message: 503
       From: DESCRIPTION
```
Please note these two DOIs are valid, both are from the journal Applied Psychological Measurement. I know it is unusual to upload a new version so fast but as I indicated above and I indicated by mail I could really use this new version for a course I am going to teach this April 10th. Sorry for the inconvenience. Finally I am using the doSNOW package, to be able to show a progress bar in the advance of the calculations

## cdcatR 1.0.4 release (March 31, 2022)
* CRAN Package Check Results found the following issue "Calling && or || with either argument of length greater than one now gives a warning (which it is intended will become an error).". This has been corrected by removing the problematic cases and contacting the author of the simGDINA function. Thanks Prof Brian Ripley and CRAN. 
* We have also included minor changes in att.plot and LR.2step functions

* cdcatR 1.0.4 re-submission:

#### Test envirnoments
* Local Windows 10 x64, R 4.1.2
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

#### R CMD check results
* There were no ERRORs or WARNINGs.
* There were 2 NOTES:
* First note: 
```
Found the following files/directories: 'lastMiKTeXException'
```
Please note that this might be due to a bug/crash in MiKTeX, as is noted in R-hub issue #503, which can be found at r-hub/rhub#503

* Second note:
```
Maintainer: 'Miguel A. Sorrel <miguel.sorrel@uam.es>'
   Found the following (possibly) invalid DOIs:
     DOI: 10.1177/0146621614554650
       From: DESCRIPTION
       Status: Service Unavailable
     DOI: 10.1177/0146621618813113
       Status: Service Unavailable
       Message: 503
       From: DESCRIPTION
```
Please note these two DOIs are valid, both are from the journal Applied Psychological Measurement

## cdcatR 1.0.3 release (July 6, 2021)
* We have included new arguments for the cdcat and gen.itembank functions of the package

## cdcatR 1.0.2 release (September 9, 2020)
* Fixed minor bugs for the cdcat() function
* We get this note: "> checking for future file timestamps ... NOTE
  unable to verify current time". I believe this is related to the check function 

## cdcatR 1.0.1 release (July 1, 2020)
* CRAN Package Check Results found an issue for noLD. We have corrected this [if (sum(probs) != 1) {] has been replaced]. Thanks Prof Brian Ripley and CRAN

## Second resubmission (June 8, 2020)
* We have addressed the only concern mentioned in the mail we got on June 2, 2020 (thanks Swetlana Herbrandt). Specifically, 
* We have includeda small executable (not wrapped in \dontrun or \donttest) example for the main function (i.e., cdcat()). It is wrapped in \dontshow so that this informal test is not shown in the help page
* We now use \donttest instead of \dontrun for examples that take longer than 10 seconds. We have checked that all the examples work fine  

## First sesubmission (May 26, 2020)
* We have addressed the two concerns mentioned in the mail we got on May 25, 2020 (thanks Jelena Saf). Specifically, 
* We have removed the redundant "in R" from the title
* A more elaborate package description was added in the description field of the DESCRIPTION file. Some references were included
* TRUE and FALSE are now always used instead of T and F. No variables are named T or F 

## Test environments
* local Windows installation, R 3.6
* Ubuntu Linux 16.04 LTS, R-release, GCC (check_rhub)
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit (check_rhub)
* Fedora Linux, R-devel, clang, gfortran (check_rhub)
* (check_win_devel)

## R CMD check results
There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* Maintainer: ‘Miguel A. Sorrel <miguel.sorrel@uam.es>’
  
  New submission
  Possibly mis-spelled words in DESCRIPTION:
    Barrada (7:277)
    Chiu (7:371)
    Abad (7:689)
    Kaplan (7:252)
    Olea (7:698)
    Tsai (7:381)
    de (7:260, 7:676)
    
* Please note thate these words are surnames. The spelling is correct
