EXE=raremetalTest
TOOLBASE = af1KG annoGroups condAnalysis Meta MetaUtility WritePDF SummaryFileReader metaElement
SRCONLY = Main.cpp 
#SRCONLY = test.cpp
HDRONLY = My_mvt.h qfc.h Rmath.h WriteLog.h imhof.h Rmath.h

LIB_RAREMETAL_PATH = ../../libRareMetal/
USER_INCLUDES = -I ../../otherLib/eigen-3.2.0 -I$(LIB_RAREMETAL_PATH)include

# Set additional values prior to including the base Makefiles.
# Set dependencies
opt: $(LIB_RAREMETAL_PATH)libRareMetal.a
debug: $(LIB_RAREMETAL_PATH)libRareMetal_debug.a
profile: $(LIB_RAREMETAL_PATH)libRareMetal_profile.a

# Define USER_LIBS as it should be used in each build type
opt: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal.a
debug: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_debug.a
profile: USER_LIBS += $(LIB_RAREMETAL_PATH)libRareMetal_profile.a

OTHER_LIBS = -lRmath

PARENT_MAKE = Makefile.src
include ../../Makefile.inc

########################
# Add any additional targets here.

$(LIB_RAREMETAL_PATH)libRareMetal.a:
	echo here
	$(MAKE) -C $(@D) --no-print-directory opt
$(LIB_RAREMETAL_PATH)libRareMetal_debug.a:
	echo here
	$(MAKE) -C $(@D) --no-print-directory debug
$(LIB_RAREMETAL_PATH)libRareMetal_profile.a:
	$(MAKE) -C $(@D) --no-print-directory profile

USER_COMPILE_VARS = -DVERSION="\"${VERSION}\"" -g -O0 -lstdc++
#USER_COMPILE_VARS = -DVERSION="\"${VERSION}\"" -g -O0 -fstack-protector-all -static
binary: USER_COMPILE_VARS +=
binary: opt
