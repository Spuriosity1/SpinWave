# Makefile to compile and link QSI simulator programs
#
# Created on 01/08/2018
# Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Creative Commons Attribution License (CC-BY),
# version 4.0 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Please find a copy of the Creative Commons CC-BY License on
# <https://creativecommons.org/licenses/by/4.0/>.

# ---------- Compiler and linker directives ----------
CXX = g++-12

CXXFLAGS = -std=c++17 -MD -g -Wall -pedantic
CXXOPTS = 
# LDFLAGS = -Wl,--gc-sections

# mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# VPATH += hartreefock numerical ode spline_basis struct_hartreefock struct_rateeq

# Output directory after linking
SRC_PATH = src
BIN_PATH = bin
BUILD_PATH = build

# ------- System Include  --------------
# Change ~/include as necessary for your architecture
INC_FLAGS = -I/opt/homebrew/include

# -------- Establishing objects and sources ------

# SOURCES = $(wildcard struct/*.cpp) calc_rates.cpp
# OBJECTS = $(addprefix $(BUILD_PATH)/,$(SOURCES:.cpp=.o))

# ------------- Specific Include -----------
SUBDIRS = struct

# Add all of the subdirectories
INC = $(addprefix -I,$(SRC_PATH)/$(SUBDIRS))

# Object dependence
json_parse_test: $(BUILD_PATH)/struct/lattice.o

# ---------- Libraries needed for linking ----------

# Change or remove ~/lib as necessary for your architecture
LIB_PATH = -L/opt/homebrew/lib -larmadillo


# Generic .cc -> .o rule
$(BUILD_PATH)/%.o: $(SRC_PATH)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(CXXOPTS) $(INC_FLAGS) $(INC) $< -c -o $@

# ---------- Linking rules ----------

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
%: $(BUILD_PATH)/%.o
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(LDFLAGS) $^ \
	$(LIB_PATH) -o $(BIN_PATH)/$@


# Don't allow direct compilation 
%: %.cpp
# Don't allow removal of object files in chained implicit linking
.SECONDARY:
# get rid of all builtin rules
.SUFFIXES:

# Cleanup:
.PHONY: clean build_directories
clean:
	-rm -r $(BUILD_PATH)/**

