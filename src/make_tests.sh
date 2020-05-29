#!/bin/bash
./make.sh test test_rectangle line_rectangle constant
mv ../bin/RunAndTumble_test_test_rectangle_line_rectangle_constant ../tests/RunAndTumble_test_test_rectangle_line_rectangle_constant
../src/make.sh test test_ring line_ring constant
mv ../bin/RunAndTumble_test_test_ring_line_ring_constant ../tests/RunAndTumble_test_test_ring_line_ring_constant
