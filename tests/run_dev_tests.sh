#!/usr/bin/env bash

# DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# src/kinfin.py -g example/OrthologousGroups.txt -c example/config.txt -s example/SequenceIDs.txt -t example/tree.nwk -o example/test -p example/SpeciesIDs.txt -a example/fasta/ -f example/functional_annotation.txt --min_proteomes 2
# #src/kinfin.py -g example/OrthologousGroups.txt -c example/config.txt -s example/SequenceIDs.txt -t example/tree.nwk -o example/test -p example/SpeciesIDs.txt -a example/fasta/ -f example/functional_annotation.txt --min_proteomes 2 --test kruskal
# #src/kinfin.py -g example/OrthologousGroups.txt -c example/config.txt -s example/SequenceIDs.txt -t example/tree.nwk -o example/test -p example/SpeciesIDs.txt -a example/fasta/ -f example/functional_annotation.txt --min_proteomes 2 --test ks

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "../$DIR"

# To exit on error
set -e 

handle_error() {
    echo "Error: $1" >&2
    exit 1
}

# Function to check if a directory exists and is not empty
function is_directory_not_empty {
    local dir="$1"
    if [ -d "$dir" ] && [ "$(ls -A $dir)" ]; then
        return 0  # Directory exists and is not empty
    else
        return 1  # Directory does not exist or is empty
    fi
}


if ! is_directory_not_empty ".test_data"; then
    echo "Extracting test data..."
    tar -xzvf ./tests/test_data.tar.gz -C "./" || handle_error "Failed to extract test data."
else
    echo "Test data is already extracted and present."
fi

echo "Running basic analysis with old tool (kinfin.py)..."
src/kinfin.py -g ".test_data/basic/input/Orthogroups.txt" -c ".test_data/basic/input/kinfin.config.basic.txt" -s ".test_data/basic/input/kinfin.SequenceIDs.txt" -o "result/basic.cli.old" || handle_error "Failed to run basic analysis with old tool."

echo "Running basic analysis with new tool (main.py)..."
src/main.py analyse -g ".test_data/basic/input/Orthogroups.txt" -c ".test_data/basic/input/kinfin.config.basic.txt" -s ".test_data/basic/input/kinfin.SequenceIDs.txt" -o "result/basic.cli.new" || handle_error "Failed to run basic analysis with new tool."

echo "Comparing output of old and new tools for basic analysis..."
pytest -v ./tests/test_output_match.py --expected result/basic.cli.old --generated result/basic.cli.new

# Check pytest exit status
if [ $? -ne 0 ]; then
    echo "Basic test failed. Stopping execution."
    exit 1
fi

# If we get here, the basic test passed, so continue with advanced analysis
echo "Running advanced analysis with old tool (kinfin.py)..."
src/kinfin.py -g ".test_data/advanced/input/Orthogroups.txt" -c ".test_data/advanced/input/kinfin.config.advanced.txt" -s ".test_data/advanced/input/kinfin.SequenceIDs.txt" -o "result/advanced.cli.old" -p ".test_data/advanced/input/kinfin.SpeciesIDs.txt" -a ".test_data/advanced/input/fastas/" -t ".test_data/advanced/input/kinfin.tree.nwk" -f ".test_data/advanced/input/kinfin.functional_annotation.txt" || handle_error "Failed to run advanced analysis with old tool."

echo "Running advanced analysis with new tool (main.py)..."
src/main.py analyse -g ".test_data/advanced/input/Orthogroups.txt" -c ".test_data/advanced/input/kinfin.config.advanced.txt" -s ".test_data/advanced/input/kinfin.SequenceIDs.txt" -o "result/advanced.cli.new" -p ".test_data/advanced/input/kinfin.SpeciesIDs.txt" -a ".test_data/advanced/input/fastas/" -t ".test_data/advanced/input/kinfin.tree.nwk" -f ".test_data/advanced/input/kinfin.functional_annotation.txt" || handle_error "Failed to run advanced analysis with new tool."

echo "Comparing output of old and new tools for advanced analysis..."
pytest -v ./tests/test_output_match.py --expected result/advanced.cli.old --generated result/advanced.cli.new