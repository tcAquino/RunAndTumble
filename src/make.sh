#!/bin/bash
sed -i '' "s/.*using namespace model_.*/	using namespace model_$1;/" "RunAndTumble.cpp"
sed -i '' "s/.*using namespace domain_.*/  using namespace domain_$2;/" "RunAndTumble.cpp"
sed -i '' "s/.*using namespace initial_condition_particles_.*/  using namespace initial_condition_particles_$3;/" "RunAndTumble.cpp"
sed -i '' "s/.*using namespace initial_condition_fields_.*/  using namespace initial_condition_fields_$4;/" "RunAndTumble.cpp"
make RunAndTumble
mv "RunAndTumble" "../bin/RunAndTumble_$1_$2_$3_$4"
