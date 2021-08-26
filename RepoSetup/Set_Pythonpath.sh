#!/usr/bin/env bash


function update_path() {
    updatedPath=${!1}
    for mrd in "${@:2}"; do
	dirToAdd="${SCRIPT_DIR}/${mrd}"
	if [[ ":${updatedPath}:" != *":${dirToAdd}:"* ]]; then
	    updatedPath="${dirToAdd}:${updatedPath}"
	    echo "Added ${dirToAdd} to $1."
	fi
    done
    [[ "${updatedPath: -1}" == ":" ]] && updatedPath="${updatedPath::-1}"
    export "${1}=${updatedPath}"
    echo "Current env variable ${1}:"
    printenv "$1"
}


# One could also loop over the env variables: declare -a PY_ENV_VARS=("PYTHONPATH" "JUPYTER_PATH") ...
declare -a MODULE_REL_DIRS=("..")

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

update_path "PYTHONPATH" "${MODULE_REL_DIRS[@]}"
update_path "JUPYTER_PATH" "${MODULE_REL_DIRS[@]}"
