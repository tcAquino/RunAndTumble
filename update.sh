#!/usr/bin/env sh

SRC=~/Dropbox/proj
PTH=./include
XCD=~/Dropbox/xcode/projects/RunAndTumble/RunAndTumble
PMN=./src

rm -rf $PTH

mkdir -p $PTH/Bacteria
cp $SRC/Bacteria/Domains.h $PTH/Bacteria/Domains.h
cp $SRC/Bacteria/InitialConditions.h $PTH/Bacteria/InitialConditions.h
cp $SRC/Bacteria/Models.h $PTH/Bacteria/Models.h
cp $SRC/Bacteria/Output.h $PTH/Bacteria/Output.h

mkdir -p $PTH/Field
cp $SRC/Field/ScalarField.h $PTH/Field/ScalarField.h
cp $SRC/Field/Measurer.h $PTH/Field/Measurer.h

mkdir -p $PTH/general
cp $SRC/general/Constants.h $PTH/general/Constants.h
cp $SRC/general/Modular.h $PTH/general/Modular.h
cp $SRC/general/MultiArray.h $PTH/general/MultiArray.h
cp $SRC/general/Operations.h $PTH/general/Operations.h
cp $SRC/general/Ranges.h $PTH/general/Ranges.h
cp $SRC/general/useful.h $PTH/general/useful.h

mkdir -p $PTH/Geometry
cp $SRC/Geometry/Shape.h $PTH/Geometry/Shape.h
cp $SRC/Geometry/Coordinates.h $PTH/Geometry/Coordinates.h

mkdir -p $PTH/Grid
cp $SRC/Grid/Grid.h $PTH/Grid/Grid.h

mkdir -p $PTH/ODE
cp $SRC/ODE/CrankNicolson.h $PTH/ODE/CrankNicolson.h
cp $SRC/ODE/BoundaryConditions.h $PTH/ODE/BoundaryConditions.h

mkdir -p $PTH/Stochastic
mkdir -p $PTH/Stochastic/CTRW
cp $SRC/Stochastic/CTRW/Boundary.h $PTH/Stochastic/CTRW/Boundary.h
cp $SRC/Stochastic/CTRW/CTRW.h $PTH/Stochastic/CTRW/CTRW.h
cp $SRC/Stochastic/CTRW/Grid.h $PTH/Stochastic/CTRW/Grid.h
cp $SRC/Stochastic/CTRW/JumpGenerator.h $PTH/Stochastic/CTRW/JumpGenerator.h
cp $SRC/Stochastic/CTRW/Measurer.h $PTH/Stochastic/CTRW/Measurer.h
cp $SRC/Stochastic/CTRW/Particle.h $PTH/Stochastic/CTRW/Particle.h
cp $SRC/Stochastic/CTRW/PTRW.h $PTH/Stochastic/CTRW/PTRW.h
cp $SRC/Stochastic/CTRW/Reaction.h $PTH/Stochastic/CTRW/Reaction.h
cp $SRC/Stochastic/CTRW/State.h $PTH/Stochastic/CTRW/State.h
cp $SRC/Stochastic/CTRW/StateGetter.h $PTH/Stochastic/CTRW/StateGetter.h
cp $SRC/Stochastic/CTRW/StateSwitcher.h $PTH/Stochastic/CTRW/StateSwitcher.h
cp $SRC/Stochastic/CTRW/TimeGenerator.h $PTH/Stochastic/CTRW/TimeGenerator.h
cp $SRC/Stochastic/CTRW/Transitions_State.h $PTH/Stochastic/CTRW/Transitions_State.h

mkdir -p $PMN
cp $XCD/main.cpp $PMN/RunAndTumble.cpp
