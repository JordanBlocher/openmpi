CC := mpic++
SRCDIR := src
INCLUDEDIR := include
UTIL := util
BUILDDIR := build
CFLAGS := -Wall -std=c++0x -g
LIBS :=  -lm
BINDIR := bin
SRCEXT := cpp

all: bin/summa bin/cannon

bin/cannon: src/cannon.cpp 
	$(CC) src/cannon.cpp -I$(INCLUDEDIR) $(CFLAGS) -I$(UTIL) $(TARGET) $(LIBS) -o bin/cannon

bin/summa: src/summa.cpp
	$(CC) src/summa.cpp -I$(INCLUDEDIR) $(CFLAGS) -I$(UTIL) $(TARGET) $(LIBS) -o bin/summa
 
clean:
	@echo " Cleaning..."; rm bin/summa bin/cannon

-include $(DEPS)

.PHONY: clean


