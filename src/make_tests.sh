#!/bin/bash
echo "Compiling rectangle test..."
./make.sh test test_rectangle line_rectangle constant
mv ../bin/RunAndTumble_test_test_rectangle_line_rectangle_constant ../tests/RunAndTumble_test_rectangle
echo "Done!"
echo "Compiling ring test..."
./make.sh test test_ring line_ring constant
mv ../bin/RunAndTumble_test_test_ring_line_ring_constant ../tests/RunAndTumble_test_ring
echo "Done!"
echo "Compiling cif test..."
./make.sh test cif_rectangle line_rectangle constant
mv ../bin/RunAndTumble_test_cif_rectangle_line_rectangle_constant ../tests/RunAndTumble_test_cif
echo "Done!"
