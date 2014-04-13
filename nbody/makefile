CC := mpic++
SRCDIR := src
INCLUDEDIR := include
UTIL := util
BUILDDIR := build
CFLAGS := -Wall -std=c++0x -g
LIBS :=  -lm
TARGET := bin/nbody
PNG := png
DATA := data
 
SRCEXT := cpp
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

DEPS := $(OBJECTS:.o=.deps)
 
$(TARGET): $(OBJECTS)
	@echo " Linking..."; $(CC) $^ -o $(TARGET) $(LIBS)
 
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT) 
	@mkdir -p $(BUILDDIR) -p $(PNG) -p $(DATA)
	@echo " CC $<"; $(CC) $(CFLAGS) -I$(INCLUDEDIR) -I$(UTIL) -MD -MF $(@:.o=.deps) -o $@ -c $< $(LIBS)
 
clean:
	@echo " Cleaning..."; $(RM) -r $(BUILDDIR) $(DATA) $(PNG) $(TARGET)

-include $(DEPS)

.PHONY: clean


