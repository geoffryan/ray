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

DEBUG    = -g -Wall
OPT      = -O3
INCLUDES = 
CFLAGS   = -c $(INCLUDES)
LDFLAGS  = -framework Accelerate -lm -lz 

ifeq ($(strip $(USE_DEBUG)), 1)
	CFLAGS += $(DEBUG)
endif
ifeq ($(strip $(USE_OPT)), 1)
	CFLAGS += $(OPT)
endif
ifeq ($(strip $(USE_OMP)), 1)
	CFLAGS += -DUSEOMP $(OMP_FLAG)
	LDFLAGS += $(OMP_FLAG)
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
