##############################################################################
# (c) The copyright relating to this work is owned jointly by the Crown,
# Met Office and NERC 2014. However, it has been created with the help of the
# GungHo Consortium, whose members are identified at
# https://puma.nerc.ac.uk/trac/GungHo/wiki
##############################################################################

include ../../Makefile.inc

BIN_DIR   = ../../bin

EXE  = dynamo
-include $(OBJ_DIR)/$(EXE).mk

.PHONY: all
all: $(BIN_DIR)/$(EXE)

$(BIN_DIR)/$(EXE): $(OBJ_DIR)/$(EXE) | $(BIN_DIR)
	@echo "Installing $@"
	@cp $< $(BIN_DIR)

# Directories
$(BIN_DIR):
	@echo "Creating $@"
	@mkdir -p $@

$(OBJ_DIR):
	@echo "creating $@"
	@mkdir -p $@

# Rules
$(OBJ_DIR)/%.mod: $(OBJ_DIR)/%.o
	@echo "Require $@"

$(OBJ_DIR)/%.o: %.F90 | $(OBJ_DIR)
	@echo "Compile $<"
	$(FC) $(FFLAGS) $(FWARNINGFLAGS) $(F_MOD_DESTINATION_ARG) \
	      -I $(OBJ_DIR) -c -o $@ $<

$(OBJ_DIR)/%.o: %.f90 | $(OBJ_DIR)
	@echo "Compile $<"
	$(FC) $(FFLAGS) $(FWARNINGFLAGS) $(F_MOD_DESTINATION_ARG) \
	      -I $(OBJ_DIR) -c -o $@ $<

$(OBJ_DIR)/$(EXE): $($(shell echo $(EXE) | tr a-z A-Z)_OBJS)
	@echo "Linking $@"
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $^

.PHONY: clean
clean:
	-rm -rf $(OBJ_DIR)

-include $(OBJ_DIR)/$(DEP_FILE)
