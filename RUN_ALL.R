# RUN_ALL.R — published-panel model verification, then appendix sensitivity checks.
# DO NOT mistake this for a verified original-source-to-FINAL30 pipeline.
PROJECT <- 'C:/Users/USERNAME/Downloads/JPR_REPLICATION_RELEASE_FROM_ACTUAL_CODE'
source(file.path(PROJECT,'INSTALL_LOCAL.R'),local=FALSE)
source(file.path(PROJECT,'RUN_PUBLICATION.R'),local=FALSE)
source(file.path(PROJECT,'RUN_APPENDIX.R'),local=FALSE)
cat('\nMODEL AND APPENDIX REPRODUCTION TIER FINISHED.\n')
cat('Full raw-to-published FINAL30 reconstruction still awaits the missing original FINAL30 construction script.\n')
