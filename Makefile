TARGET = main
CPP = g++
CPPFLAGS = -c -std=c++11 -Wall -Wextra -Wpedantic -O3 -MD

SRCPATH = ./src
OBJDIR_RELEASE = obj/Release
OBJ_RELEASE = \
	$(OBJDIR_RELEASE)/main.o \
	$(OBJDIR_RELEASE)/timer.o \
	$(OBJDIR_RELEASE)/test.o 

HEADERS = \
	$(SRCPATH)/trie_based.h \
	$(SRCPATH)/timer.h \
	$(SRCPATH)/data_io.h \
	$(SRCPATH)/quantile.h \
	$(SRCPATH)/test.h 

all: release

clean: clean_release

release: before_release out_release

before_release:
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

out_release: $(OBJ_RELEASE) $(HEADERS)
	$(CPP) -o $(TARGET) $(OBJDIR_RELEASE)/*.o

$(OBJDIR_RELEASE)/main.o: $(SRCPATH)/main.cpp
	$(CPP) $(CPPFLAGS) $(SRCPATH)/main.cpp -o $(OBJDIR_RELEASE)/main.o
$(OBJDIR_RELEASE)/timer.o: $(SRCPATH)/timer.cpp
	$(CPP) $(CPPFLAGS) $(SRCPATH)/timer.cpp -o $(OBJDIR_RELEASE)/timer.o
$(OBJDIR_RELEASE)/test.o: $(SRCPATH)/test.cpp
	$(CPP) $(CPPFLAGS) $(SRCPATH)/test.cpp -o $(OBJDIR_RELEASE)/test.o

clean_release:
	rm $(OBJDIR_RELEASE)/*.o
	rm $(OBJDIR_RELEASE)/*.d
	rm $(TARGET)

-include $(OBJ_RELEASE:.o=.d)