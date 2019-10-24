#!/usr/bin/env bash
ln -sf ../../GPR/GPR/lib/ .
export PROJECT_LFLAGS="$PROJECT_LFLAGS -L./lib -ltecio"
