CXX      = g++
CXXFLAGS = -std=c++17 -O3 -march=native -Wall -Wextra -DNDEBUG
LDFLAGS  = -lm

SRC_DIR  = src
OBJ_DIR  = obj
TARGET   = aero3d

SOURCES  = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS  = $(SOURCES:$(SRC_DIR)/%.cpp=$(OBJ_DIR)/%.o)

.PHONY: all clean debug

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

debug: CXXFLAGS = -std=c++17 -O0 -g -Wall -Wextra
debug: clean $(TARGET)

clean:
	rm -rf $(OBJ_DIR) $(TARGET)
