#!/bin/sh
# script to format all sources
# in subdirectories according 
# to .clang-format - 2019 - Haroldo
# modified from Cbc - Simon

echo formatting all source using .clang-format
find ./ -iname '*.h' -o -iname '*.cpp' -exec clang-format -style=llvm -i {} +

# Check for formatting. Use this in CI pipeline
# find ./ -iname '*.h' -o -iname '*.cpp' -exec clang-format -style=llvm -output-replacements-xml {} + | grep -c "<replacement " | grep 0 >/dev/null
