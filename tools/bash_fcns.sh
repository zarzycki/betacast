#!/bin/bash

# Strip surrounding quotes from string [$1: variable name]
function strip_quotes() {
  local -n var="$1"
  [[ "${var}" == \"*\" || "${var}" == \'*\' ]] && var="${var:1:-1}"
}
export -f strip_quotes
