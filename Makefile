### Makefile for compilation of analysis macros
### to be copied into the analysis directory or used with make -f <Makefile location>
### and perhaps edited to correct the FWK_BASE_DIR e.g. make example_gen_ttbar FWK_BASE_DIR=<fwk base dir> -f <Makefile location>

FWK_BASE_DIR = ..
COMPILER = g++
SAVE_TEMPS = false
DEBUG = false
FWK_SRC_DIR = $(FWK_BASE_DIR)/src/
FWK_PLUGINS_DIR = $(FWK_BASE_DIR)/plugins/

CXXFLAGS = $(shell root-config --cflags --evelibs) -std=c++17 -O3 -Wall -Wextra -Wpedantic -Werror -Wno-float-equal -Wno-sign-compare -Wfatal-errors
FWK_INCLUDES = -I $(FWK_SRC_DIR) -I $(FWK_PLUGINS_DIR)


%: %.cc
	$(eval CXXFLAGS += $(if $(filter $(COMPILER), clang++), -stdlib=libstdc++ -frelaxed-template-template-args, ))
	$(eval CXXFLAGS += $(if $(filter $(SAVE_TEMPS), true), -save-temps, ))
	$(eval CXXFLAGS += $(if $(filter $(DEBUG), true), -ggdb3, ))
	$(info $(COMPILER) $(CXXFLAGS) $(FWK_INCLUDES) $< -o $@)
	@ $(COMPILER) $(CXXFLAGS) $(FWK_INCLUDES) $< -o $@
