
commit 8eea110ffd1c2b71f0d0fd58fb144baf08ea3002 (HEAD -> master, origin/master)
Author: Blahak Ulrich <ulrich.blahak@dwd.de>
Date:   Wed Apr 8 13:32:05 2020 +0200

-  Initial creation of this repository.

-  It contains the subset of sources from the DACE code which are necessary
   for the data assimilation online forward operators in ICON.

-  ICON uses this repository and its master branch as a submodule "dace_icon".

-  ICON expects these DACE files under the subdirectory ./src_for_icon/.
   From ICON side, this directory would in principle be allowed to have subdirectories,
   but this is not the case at the moment.

