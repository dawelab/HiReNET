# =========================================================
# HiReNET: Installation and Utility Commands
# =========================================================

# Allow users to pick a prefix: e.g. `make install PREFIX=/opt/hirenet`
PREFIX ?= $(HOME)/HiReNET
BIN_DIR := $(PREFIX)/bin
SCRIPT_DIR := $(PREFIX)/scripts
R_DIR := $(PREFIX)/R

.PHONY: all install dispatcher help clean perms

all: help

# Copy the repo tree into $(PREFIX) if needed, then build dispatcher there
install:
	@mkdir -p "$(BIN_DIR)" "$(SCRIPT_DIR)" "$(R_DIR)"
	# Copy scripts/R if this isn't already the live tree
	@if [ ! -f "$(SCRIPT_DIR)/00_runALL.sh" ]; then \
	  echo " Copying scripts/ to $(SCRIPT_DIR)"; \
	  cp -r scripts/* "$(SCRIPT_DIR)/"; \
	fi
	@if [ ! -f "$(R_DIR)/S1_network_unmerge_bin_noplot.R" ]; then \
	  echo " Copying R/ to $(R_DIR)"; \
	  cp -r R/* "$(R_DIR)/"; \
	fi
	# Place/refresh dispatcher under $(BIN_DIR)
	@echo " Building HiReNET command dispatcher in $(BIN_DIR)..."
	@ROOT="$(PWD)"; \
	  sed 's|^ROOT=.*|ROOT="$${ROOT}"|' make_dispatcher.sh >/dev/null 2>&1 || true; \
	  bash make_dispatcher.sh
	@chmod +x "$(BIN_DIR)/HiReNET"
	@$(MAKE) perms
	@echo ""
	@echo " HiReNET installed successfully to $(PREFIX)"
	@echo "Add to PATH (if not already):"
	@echo "    echo 'export PATH=\"$(BIN_DIR):\$$PATH\"' >> ~/.bashrc && source ~/.bashrc"
	@echo "    # or zsh: echo 'export PATH=\"$(BIN_DIR):\$$PATH\"' >> ~/.zshrc && source ~/.zshrc"
	@echo "You can now run: HiReNET --help"

# Rebuild dispatcher only (uses repo-local paths)
dispatcher:
	@echo " Building HiReNET command dispatcher..."
	bash make_dispatcher.sh
	@$(MAKE) perms
	@echo ""

# Make scripts executable safely (no-glob failures)
perms:
	@echo " Ensuring scripts are executable..."
	@find "$(SCRIPT_DIR)" -type f -name '*.sh' -exec chmod +x {} \; 2>/dev/null || true
	@echo " Made scripts executable under $(SCRIPT_DIR)"

help:
	@echo "Usage:"
	@echo "  make install         Copy (if needed), build dispatcher, fix perms"
	@echo "  make dispatcher      Rebuild dispatcher only"
	@echo "  make clean           Remove dispatcher"

clean:
	@echo " Cleaning up..."
	@rm -f "$(BIN_DIR)/HiReNET" || true
	@echo "Removed $(BIN_DIR)/HiReNET"