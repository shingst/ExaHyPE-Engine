## Contributing to mexa

The `mexa` code is currently maintained withing the 
[ExaHyPE-Engine](https://gitlab.lrz.de/exahype/ExaHyPE-Engine) repository.
However, a public
[standalone version](https://bitbucket.org/svek/mexa) is also provided.

### Dependencies to ExaHyPE

Despite the naming may suggest it: There are no dependencies at all from the
mexa code to ExaHyPE. The other way around it may true for certain ExaHyPE
applications: They probably include the `mexa-cpp` code and thus are dependent.

Since the `mexa-cpp` is an (almost) single-file code, it is suitable to copy
instead of linking the codes and keep them in sync.

### Merge the two repositories

Here are some commands to merge the repositories:

* Inside a **copy** of the Engine repository: `git filter-branch --prune-empty --subdirectory-filter  Miscellaneous/MetaSpecfile/  master`
* Then register the repo as a remote (branch) in the standalone repo
* Then just `git merge -X theirs exahype-engine-repository/master`
* And if uploading fails, do `git push -f origin master`

*Note*: This attemp will treat the standalone repository as recieving slave
and overwrite changes done their.

