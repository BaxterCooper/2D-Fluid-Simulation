# 2D-Fluid-Simulation

A 2D fluid simulation inspired by Jos Stam's paper [Real-Time Fluid Dynamics for Games](http://graphics.cs.cmu.edu/nsp/course/15-464/Fall09/papers/StamFluidforGames.pdf) and Mike Ash's thesis [Fluid Simulation for Dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html).

For an inspirational video on the physics and concepts used in this simulation please watch Gonkee's video [But How DO Fluid Simulations Work?](https://www.youtube.com/watch?v=qsYE1wMEMPA).

![](2D-Fluid-Simulation-Demo.gif)

## Prerequisites 

### Install SFML
```bash
sudo apt-get install libsfml-dev
```

## Compilation

### Windows
```bash
make
```

### Linux
```bash
g++ -c main.cpp && g++ main.o -o main -lsfml-graphics -lsfml-window -lsfml-system
```

## Execution
```bash
./main
```
