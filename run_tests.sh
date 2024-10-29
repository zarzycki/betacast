#!/bin/bash
# Betacast Test Runner
#
# This script runs comparison tests between Python and NCL implementations
# of the BetaCast weather model initialization system. Each test will:
#   1. Run the Python version
#   2. Run the NCL version
#   3. Compare the results
#   4. Clean up temporary files

# Colors for output
GREEN='\033[0;32m'
RED='\033[0;31m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Help function
show_help() {
    cat << EOF
Usage:
  ./run_tests.sh [options]

Options:
  --test-dir DIR    Directory containing test files (default: ./tests)
  --test NAME       Run a specific test by name (without .sh extension)
  --help           Show this help message

Example usage:
  ./run_tests.sh                         # Run all tests in ./tests
  ./run_tests.sh --test se_basic         # Run only the SE basic test
  ./run_tests.sh --test-dir /path/tests  # Run tests from specific directory
EOF
    exit 0
}

# Test results tracking
PASSED=0
FAILED=0
TOTAL=0
declare -a PASSED_TESTS=()
declare -a FAILED_TESTS=()
declare -A TEST_FAILURE_REASONS=()

# Function to run a test
run_test() {
    local test_name=$1
    local test_file=$2

    echo -e "\n${BLUE}=====================================${NC}"
    echo -e "${BLUE}Running test: ${test_name}${NC}"
    echo -e "${BLUE}=====================================${NC}\n"

    # Source the test file to get the environment variables and functions
    source "$test_file"

    # Run the Python version
    echo -e "${BLUE}Running Python version...${NC}"
    if ! run_python_test; then
        TEST_FAILURE_REASONS[$test_name]="Python test failed"
        return 1
    fi
    echo -e "${GREEN}Python test completed${NC}"

    # Run the NCL version
    echo -e "\n${BLUE}Running NCL version...${NC}"
    run_ncl_test
    ncl_status=$?
    if [ "$ncl_status" -ne 9 ]; then
        TEST_FAILURE_REASONS[$test_name]="NCL test failed (exit code: $ncl_status)"
        return 1
    fi
    echo -e "${GREEN}NCL test completed${NC}"

    # Run validation
    echo -e "\n${BLUE}Validating results...${NC}"
    if ! run_validation; then
        TEST_FAILURE_REASONS[$test_name]="Validation failed"
        return 1
    fi
    echo -e "${GREEN}Validation passed${NC}"

    # Run cleanup if defined
    if type cleanup >/dev/null 2>&1; then
        echo -e "\n${BLUE}Running cleanup...${NC}"
        cleanup
    fi

    return 0
}

# Main function
main() {
    # Default test directory
    TEST_DIR="tests"

    # Parse command line arguments
    while [[ $# -gt 0 ]]; do
        case $1 in
            --help|-h)
                show_help
                ;;
            --test-dir)
                TEST_DIR="$2"
                shift 2
                ;;
            --test)
                SPECIFIC_TEST="$2"
                shift 2
                ;;
            *)
                echo "Unknown option: $1"
                echo "Use --help for usage information"
                exit 1
                ;;
        esac
    done

    # Check if test directory exists
    if [ ! -d "$TEST_DIR" ]; then
        echo "Test directory $TEST_DIR does not exist!"
        exit 1
    fi

    # Run tests
    if [ -n "$SPECIFIC_TEST" ]; then
        # Run specific test
        if [ -f "$TEST_DIR/$SPECIFIC_TEST.sh" ]; then
            if run_test "$SPECIFIC_TEST" "$TEST_DIR/$SPECIFIC_TEST.sh"; then
                PASSED_TESTS+=("$SPECIFIC_TEST")
                ((PASSED++))
            else
                FAILED_TESTS+=("$SPECIFIC_TEST")
                ((FAILED++))
            fi
            ((TOTAL++))
        else
            echo "Test $SPECIFIC_TEST not found!"
            exit 1
        fi
    else
        # Run all tests
        for test_file in "$TEST_DIR"/*.sh; do
            if [ -f "$test_file" ]; then
                test_name=$(basename "$test_file" .sh)
                if run_test "$test_name" "$test_file"; then
                    PASSED_TESTS+=("$test_name")
                    ((PASSED++))
                else
                    FAILED_TESTS+=("$test_name")
                    ((FAILED++))
                fi
                ((TOTAL++))
            fi
        done
    fi

    # Print summary
    echo -e "\n${BLUE}=====================================${NC}"
    echo -e "${BLUE}Test Summary${NC}"
    echo -e "${BLUE}=====================================${NC}"
    echo -e "Total tests: $TOTAL"
    echo -e "Passed: ${GREEN}$PASSED${NC}"
    echo -e "Failed: ${RED}$FAILED${NC}"

    # Print passed tests
    if [ ${#PASSED_TESTS[@]} -gt 0 ]; then
        echo -e "\n${GREEN}Passed Tests:${NC}"
        for test in "${PASSED_TESTS[@]}"; do
            echo -e "${GREEN}✓ $test${NC}"
        done
    fi

    # Print failed tests with reasons
    if [ ${#FAILED_TESTS[@]} -gt 0 ]; then
        echo -e "\n${RED}Failed Tests:${NC}"
        for test in "${FAILED_TESTS[@]}"; do
            echo -e "${RED}✗ $test${NC}"
            echo -e "  Reason: ${TEST_FAILURE_REASONS[$test]}"
        done
    fi

    # Return non-zero if any tests failed
    [ "$FAILED" -eq 0 ]
}

# Run main function
main "$@"