TARGET	= exe
CC		= g++
RM		= rm -f
CFLAGS	= -g -Wall
LIBS	= -lgsl -lgslcblas

default: $(TARGET)
all: default

OBJECTS = $(patsubst %.cpp, %.o, $(wildcard *.cpp))
HEADERS = $(wildcard *.h)

%.o: %.cpp $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

.PRECIOUS: $(TARGET) $(OBJECTS)

$(TARGET): $(OBJECTS)
	$(CC) $(OBJECTS) -Werror -Wall $(LIBS) -o $@

clean:
	$(RM) *.o
	$(RM) $(TARGET)
