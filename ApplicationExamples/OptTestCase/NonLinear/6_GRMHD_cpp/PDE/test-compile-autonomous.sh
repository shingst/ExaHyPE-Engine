#!/bin/bash

mv test.cpp{-disable,}
g++ -otest_autonomous --std=c++11 -Wall -DTEST_NEW_PDE_AUTONOMOUSLY test.cpp
mv test.cpp{,-disable}
