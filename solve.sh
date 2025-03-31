#!/bin/sh

MPSFILE=$1
SOLFILE=${MPSFILE%.[Mm][Pp][Ss]}.sol

LSMIP_PATH=`dirname "$0"`
CONVERT_SCRIPT_PATH=${LSMIP_PATH}/convert.py
TEMP_PATH=${LSMIP_PATH}/tmp/
echo Solving "$MPSFILE" "->" "$SOLFILE"
${LSMIP_PATH}/lsmip "$MPSFILE" "$SOLFILE" --convert-script-path "$CONVERT_SCRIPT_PATH" --temp-path "$TEMP_PATH"
