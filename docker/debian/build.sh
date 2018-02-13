#!/usr/bin/env bash

set -x

docker build --no-cache -t debian-esl:full  . 
