#!/bin/bash

# Usage:
#   1. Run in the current directory:
#      ./get-py-imports.sh
#
#   2. Specify a directory to search:
#      ./get-py-imports.sh /path/to/directory
#
#   3. Search multiple directories:
#      ./get-py-imports.sh /path/to/dir1 /path/to/dir2
#
#   4. Use with other commands:
#      ./get-py-imports.sh | wc -l  # Count unique imports
#      ./get-py-imports.sh | grep numpy  # Check for numpy imports
#
# Note: Ensure the script has execute permissions (chmod +x get-py-imports.sh)

# Define the directory to search (use the current directory by default)
SEARCH_DIR=${1:-.}

### CMZ: OLD CODE
# # Find all Python files in the directory
# find "$SEARCH_DIR" -name "*.py" | while read -r file; do
#   echo "Processing $file"
#   # Extract all the imported libraries from each Python file
#   grep -oP '^\s*(import|from)\s+\K[\w.]+' "$file"
# done | sort | uniq

# Find all Python files in the directory and extract unique imports
find "$SEARCH_DIR" -name "*.py" -exec grep -hoP '^\s*(import|from)\s+\K[\w.]+' {} + | sort -u