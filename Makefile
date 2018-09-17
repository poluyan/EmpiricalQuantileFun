TARGET = main
CPP = g++
CPPFLAGS = -c -std=c++11 -Wall -Wextra -O3 -MD

SRCPATH = ./src
OBJDIR_RELEASE = obj/Release
OBJ_RELEASE = \
	$(OBJDIR_RELEASE)/main.o 


all: release

clean: clean_release

release: before_release out_release

before_release:
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

out_release: $(OBJ_RELEASE) 
	$(CPP) -o $(TARGET) $(OBJDIR_RELEASE)/*.o

$(OBJDIR_RELEASE)/main.o: $(SRCPATH)/main.cpp
	$(CPP) $(CPPFLAGS) $(SRCPATH)/main.cpp -o $(OBJDIR_RELEASE)/main.o


clean_release:
	rm $(OBJDIR_RELEASE)/*.o
	rm $(OBJDIR_RELEASE)/*.d
	rm $(TARGET)
