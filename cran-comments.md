## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Review 2026-02-20

Reviewer:  konstanze.lauseker@wu.ac.at

* If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file in the form

Paper in review process, will add as citationa and also in DESCRIPTION once completed. 

* Please write TRUE and FALSE instead of T and F. Please don't use "T" or "F" as vector names.

Changed all cases where T was used to TRUE in callKaryotypes() and also in ddirichletultinomial()

* You write information messages to the console that cannot be easily
suppressed.

All uses of cat have been changed to stop(), warning() or message() depending on context. 

