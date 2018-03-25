#!/bin/bash

echo "Copying makefile to tests directory..."
cp ./makeFunctionSystemTest ../tests

echo "Entering test directory..."
cd ../tests

echo "Building test..."
make -f ./makeFunctionSystemTest

echo "Cleaning up..."
rm ./makeFunctionSystemTest
rm ./*.o
