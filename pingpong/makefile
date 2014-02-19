CC := mpic++
SRCDIR := src
INCLUDEDIR := include
UTIL := util
BUILDDIR := build
CFLAGS := -Wall -g
LIBS :=  -lm
TARGET := bin/pingpong
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

DEPS := $(OBJECTS:.o=.deps)
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $^ -o $(TARGET) $(LIBS)
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(CFLAGS) -I$(INCLUDEDIR) -I$(UTIL) -MD -MF $(@:.o=.deps) -o $@ -c $< $(LIBS)
 
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(TARGET)

-include $(DEPS)

.PHONY: clean


