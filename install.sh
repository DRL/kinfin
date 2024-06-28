#!/usr/bin/env bash

# logging function
log() {
    local GREEN='\033[0;32m'
    local YELLOW='\033[0;33m'
    local RED='\033[0;31m'
    local NO_COLOR='\033[0m'

    local level=$1
    local message=$2

    case $level in
        INFO)
            echo -e "[${NO_COLOR}INFO${NO_COLOR}] - $message"
            ;;
        SUCCESS)
            echo -e "[${GREEN}SUCCESS${NO_COLOR}] - $message"
            ;;
        ERROR)
            echo -e "[${RED}ERROR${NO_COLOR}] - $message" >&2
            ;;
        *)
            echo "Invalid log level: $level"
            ;;
    esac
}

# Check dependencies exist
check_dependencies() {
    log INFO "Checking dependencies..."

    local dependencies=("wget" "gunzip")
    local missing_dependencies=()

    for dependency in "${dependencies[@]}"; do
        local item=$(command -v "$dependency")
        if [ ! -x "$item" ]; then
            missing_dependencies+=("$dependency")
        fi
    done

    if [ ${#missing_dependencies[@]} -gt 0 ]; then
        log ERROR "Missing dependencies: ${missing_dependencies[*]}. Please install them."
        exit 1
    else
        for dependency in "${dependencies[@]}"; do
            log SUCCESS "$dependency is installed."
        done
        log SUCCESS "All dependencies are installed."
        return 0
    fi
}

# Function to download a file
download_file() {
    local url=$1
    local filename=$2

    log INFO "Downloading $filename from $url"
    $(which wget) -np -nd -qN --show-progress "$url" -P "$DIR/data/"

    if [ $? -eq 0 ]; then
        log SUCCESS "Downloaded $filename"
    else
        log ERROR "Failed to download $filename from $url"
        exit 1
    fi
}

# Extract .gz files
extract_gzip() {
    local gz_file=$1
    local dest=$2
    
    log INFO "Extracting $gz_file..."

    $(which gunzip) -c "$gz_file" > "$dest"
    
    if [ $? -eq 0 ]; then
        log SUCCESS "Extracted $gz_file at $dest"
    else
        log ERROR "Failed to extract $gz_file. Please download kinfin again."
        exit 1
    fi
}



main() {
    # Set working directory
    DIR="$(cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

    check_dependencies

    log INFO "Checking input data files..."

    local pfam_dest="$DIR/data/Pfam-A.clans.tsv.gz"
    local ipr_dest="$DIR/data/entry.list"
    local go_dest="$DIR/data/interpro2go"
    local nodesdbgz="$DIR/data/nodesdb.gz"
    local nodesdb="$DIR/data/nodesdb.txt"

    if [ ! -f "$nodesdb" ]; then
        if [ -f "$nodesdbgz" ]; then
            extract_gzip "$nodesdbgz" "$nodesdb"
        else
            log ERROR "$nodesdbgz not found. Please download kinfin again."
            exit 1
        fi
    else
        log SUCCESS "$nodesdb is already present."
    fi

    if [ ! -f "$pfam_dest" ]; then
        download_file "ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz" "Pfam-A.clans.tsv.gz"
    else
        log SUCCESS "Pfam-A.clans.tsv.gz is already present."
    fi

    if [ ! -f "$ipr_dest" ]; then
        download_file "ftp.ebi.ac.uk/pub/databases/interpro/current_release/entry.list" "entry.list"
    else
        log SUCCESS "entry.list is already present."
    fi

    if [ ! -f "$go_dest" ]; then
        download_file "ftp.ebi.ac.uk/pub/databases/interpro/current_release/interpro2go" "interpro2go"
    else
        log SUCCESS "interpro2go is already present."
    fi

    log SUCCESS "All required files downloaded."

    # Create executable
    log INFO "Creating executable..."
    echo -e '#!/usr/bin/env bash\nDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"\n$DIR/src/kinfin.py "$@"' > $DIR/kinfin && chmod +x $DIR/kinfin

    # Done
    log SUCCESS "Kinfin was installed. Please run ./kinfin --help"
}

main