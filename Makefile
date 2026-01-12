# Top-level Makefile in parent directory

CXX := g++

# Common flags
CXXFLAGS_COMMON := -std=c++17 -Wall
LDFLAGS :=

# Per-phase flags
CXXFLAGS_P1 := $(CXXFLAGS_COMMON) -O3
CXXFLAGS_P2 := $(CXXFLAGS_COMMON) -O3
CXXFLAGS_P3 := $(CXXFLAGS_COMMON) -O2 -Wextra -pedantic

# Directories
PHASE1_DIR := Phase-1
PHASE2_DIR := Phase-2
PHASE3_DIR := Phase-3

# Executables
PHASE1_BIN := phase1
PHASE2_BIN := phase2
PHASE3_BIN := phase3

.PHONY: all clean

all: $(PHASE1_BIN) $(PHASE2_BIN) $(PHASE3_BIN)

# -------- Phase 1 --------
$(PHASE1_BIN):
	$(CXX) $(CXXFLAGS_P1) -I$(PHASE1_DIR) \
		-o $(PHASE1_BIN) \
		$(PHASE1_DIR)/main.cpp \
		$(PHASE1_DIR)/functions.cpp \
		$(PHASE1_DIR)/graph_utils.cpp

# -------- Phase 2 --------
$(PHASE2_BIN):
	$(CXX) $(CXXFLAGS_P2) -I$(PHASE2_DIR) \
		-o $(PHASE2_BIN) \
		$(PHASE2_DIR)/main.cpp \
		$(PHASE2_DIR)/functions.cpp \
		$(PHASE2_DIR)/graph_utils.cpp

# -------- Phase 3 --------
$(PHASE3_BIN):
	$(CXX) $(CXXFLAGS_P3) -I$(PHASE3_DIR) \
		-o $(PHASE3_BIN) \
		$(PHASE3_DIR)/main.cpp \
		$(PHASE3_DIR)/graph.cpp \
		$(PHASE3_DIR)/scheduler.cpp

clean:
	rm -f $(PHASE1_BIN) $(PHASE2_BIN) $(PHASE3_BIN)
