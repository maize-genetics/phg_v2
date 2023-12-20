#!/bin/bash
TILEDB_URI="$HOME/temp/phgv2Tests/tempDir/testTileDBURI/"
# This adds the TILEDB_URI line to the beginning of the application.conf file
echo -e "TILEDB_URI=$TILEDB_URI\n$(cat src/main/resources/application.conf)" > src/main/resources/application.conf