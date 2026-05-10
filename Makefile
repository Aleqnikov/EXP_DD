TARGET = dd_calc

CXX = g++
CXXFLAGS = -std=c++17 -O3 -fno-fast-math -Wall -Wextra

SRCS = main.cpp DD.cpp
OBJS = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $(TARGET)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean