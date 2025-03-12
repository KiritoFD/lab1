#!/bin/bash

echo "=== Build System Check ==="

# Check for common compilation issues
echo "Checking for common compilation issues..."

# Verify include paths are correct
echo "Verifying include paths in source files..."
include_issues=0
for src_file in src/*.c; do
    if ! grep -q "../include/" "$src_file"; then
        echo "  - WARNING: $src_file may have incorrect include paths"
        include_issues=$((include_issues + 1))
    fi
done

if [ $include_issues -eq 0 ]; then
    echo "  ✓ Include paths look correct"
else 
    echo "  ! Found $include_issues files with potential include path issues"
fi

# Check for multiple main functions
echo "Checking for multiple main functions..."
main_count=$(grep -l "int main" src/*.c | wc -l)

if [ $main_count -gt 1 ]; then
    echo "  ! WARNING: Found more than one main function ($main_count found)"
    grep -l "int main" src/*.c
else
    echo "  ✓ Only one main function found"
fi

# Check if all headers in include/ are being used
echo "Verifying header usage..."
unused_headers=0
for header in include/*.h; do
    header_name=$(basename "$header")
    if ! grep -q "$header_name" src/*.[ch]; then
        echo "  - NOTE: $header_name might be unused"
        unused_headers=$((unused_headers + 1))
    fi
done

if [ $unused_headers -eq 0 ]; then
    echo "  ✓ All headers appear to be used"
fi

# Check Makefile
echo "Verifying Makefile..."
if [ -f "Makefile" ]; then
    if grep -q "SRC_DIR/core/" Makefile; then
        echo "  ! WARNING: Makefile references non-existent directory structure (SRC_DIR/core/)"
    else
        echo "  ✓ Makefile structure looks correct"
    fi
else
    echo "  ! ERROR: No Makefile found"
fi

echo ""
echo "Run 'make clean && make' to attempt compilation"
echo "=== Done ==="
