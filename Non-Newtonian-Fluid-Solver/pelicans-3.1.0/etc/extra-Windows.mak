# Common Makefile for Windows plateform 
# Settings for Extra Components

ifeq ($(LINK_PEL),1)

##### No external API
WITH_MPI      = 0
WITH_PETSc    = 0
WITH_AZTEC    = 0
WITH_SPARSKIT = 0
WITH_SILO     = 0
WITH_OPENGL   = 0
WITH_UMFPACK  = 0
WITH_METIS    = 0

endif

WITH_ZLIB = 0
WITH_EXTENDED_MATH = 1
WITH_JNI = 0


###############################################################################
ifeq ($(WITH_EXTENDED_MATH),1)
CPPFLAGS += -DEXTENDED_MATH
endif

