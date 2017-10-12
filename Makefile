# Makefile for GeoByGeo

MAKEFILE_IN = $(PWD)/Makefile.in
include $(MAKEFILE_IN)

APP      = ray

SRCEXT   = c
SRCDIR   = src
OBJDIR   = obj
BINDIR   = bin

SRCS    := $(shell find $(SRCDIR) -name '*.$(SRCEXT)')
SRCDIRS := $(shell find . -name '*.$(SRCEXT)' -exec dirname {} \; | uniq)
OBJS    := $(patsubst %.$(SRCEXT),$(OBJDIR)/%.o,$(SRCS))

PLATFORM := $(shell uname -s)

DEBUG    = -g -Wall
OPT      = -O3
INC      = -I$(H5I)
LIB      = -L$(H5L) -lhdf5
DEF      =

ifeq ($(PLATFORM),Darwin)
	DEF += -DOSX
	LIB += -framework Accelerate
endif
ifeq ($(PLATFORM),Linux)
	DEF += -DLINUX
ifeq ($(strip $(USE_MKL)), 1)
	DEF += -DMKL
	INC += -I$(MKLDIR)/include
	LIB += -L$(MKLDIR)/lib/em64t -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_lapack
else
	INC += -I$(LPI)
	LIB += -L$(LPL) -llapack
endif
endif

ifeq ($(strip $(USE_OMP)), 1)
	$(DEF) += -DUSEOMP $(OMP_FLAG)
	$(LIB) += $(OMP_FLAG)
endif

CFLAGS   = -c $(DEF) $(INC)
LDFLAGS  = $(LIB) -lm

ifeq ($(strip $(USE_DEBUG)), 1)
	CFLAGS += $(DEBUG)
endif
ifeq ($(strip $(USE_OPT)), 1)
	CFLAGS += $(OPT)
endif


.PHONY: all clean distclean


#all: $(BINDIR)/$(APP) $(VISDIR)/$(BINDIR)/$(VISAPP)
all: $(BINDIR)/$(APP)

$(VISDIR)/$(BINDIR)/$(VISAPP): .FORCE
	$(MAKE) -C $(VISDIR)

$(BINDIR)/$(APP): buildrepo $(OBJS)
	@mkdir -p `dirname $@`
	@echo "Linking $@..."
	@$(CC) $(OBJS) $(LDFLAGS) -o $@

$(OBJDIR)/%.o: %.$(SRCEXT)
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) $< -o $@

clean:
	$(RM) -r $(OBJDIR)

distclean: clean
	$(RM) -r $(BINDIR)

buildrepo:
	@$(call make-repo)

define make-repo
   for dir in $(SRCDIRS); \
   do \
	mkdir -p $(OBJDIR)/$$dir; \
   done
endef

.FORCE: 
