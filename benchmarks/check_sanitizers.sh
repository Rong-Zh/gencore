#!/usr/bin/env bash
set -eu
cd "$(dirname "$0")/.."
# Separate objects keep the normal optimized binary intact. -B ensures stale
# non-instrumented objects cannot silently bypass instrumentation.
make -B -j2 DIR_OBJ=./obj/sanitize DIR_BIN=./bin/sanitize \
    CXXFLAGS='-std=c++20 -O1 -g -MMD -MP -fsanitize=address,undefined -fno-omit-frame-pointer' \
    LDFLAGS='-fsanitize=address,undefined' ./bin/sanitize/gencore_tests
./bin/sanitize/gencore_tests
