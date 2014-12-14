#!/bin/bash

DIR="$( cd "$( dirname "$0" )" && pwd )"
cd $DIR

xcodebuild -project phycas.xcodeproj -target "Everything" -showBuildSettings > phycas_build_settings.txt
