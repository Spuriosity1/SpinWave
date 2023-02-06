#!/bin/bash

NUM_FAILED=0
NUM_PASSED=0

function check_retval {
    if [ $? ]; then
        echo "$(tput bold)$(tput setaf 2)Passed$(tput sgr0)"
        ((NUM_PASSED++))
    else
        echo "$(tput bold)$(setaf 1)FAILED$(tput sgr0)"; 
        ((NUM_FAILED++))
    fi;

}

# cd to dir enclosing script directory (i.e. c++)
cd "${0%/*}/../"
echo Building unit tests...
make json_parse_test

echo Done.
echo

# Clean up the output directory
rm test/out/*

echo -n "Json parsing... "
bin/json_parse_test test/data/square_j1j2.jsonc test/out/square_j1j2_test.jsonc > test/out/json_parse_test.log 2> test/out/json_parse_test.error
check_retval



echo
echo
echo "$(tput bold)Tests complete. ${NUM_PASSED} succeeded, ${NUM_FAILED} failed."
echo "$(tput sgr0)"
