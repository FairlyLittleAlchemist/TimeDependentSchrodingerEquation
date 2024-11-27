# Double-Slit Wave Simulation

This project simulates a quantum wave propagating through a double-slit apparatus, inspired by the famous double-slit experiment. It visualizes the evolution of the wave function over time using numerical methods and sparse matrix computations. The simulation leverages Python libraries such as `numpy`, `scipy`, and `matplotlib`.

## Features

- **Quantum Wave Function Evolution**: Models the time evolution of a quantum wave function in a 2D spatial domain.
- **Sparse Matrix Solver**: Implements the Conjugate Gradient method for efficiently solving sparse linear systems.
- **Double-Slit Potential**: Incorporates a custom double-slit potential barrier within the simulation.
- **Interactive Visualization**: Animates the modulus of the wave function to demonstrate quantum interference patterns.
- **Customizable Parameters**: Allows for adjustment of grid resolution, potential strength, and initial wave parameters.

## Simulation Overview

The simulation involves:

1. **Wave Function Initialization**: The initial wave function is a Gaussian distribution modulated by a plane wave.
2. **Sparse Matrix Construction**: Constructs sparse matrices to represent the evolution operator for time-stepping.
3. **Time Evolution**: Propagates the wave function using a Crank-Nicholson-like scheme, solved iteratively.
4. **Double-Slit Barrier**: Defines a slit-like potential that mimics the experimental setup.
5. **Visualization**: Animates the modulus of the wave function to show quantum interference.

## Code Breakdown

### Key Components

- **`psi0`**: Initializes the wave function as a Gaussian modulated by a plane wave.
- **`conjugate_gradient_sparse`**: Solves the sparse linear system using the Conjugate Gradient method.
- **Double-Slit Potential**: Creates a rectangular barrier with two slits to model the experimental setup.
- **Wave Function Modulus**: Computes the modulus (absolute value) of the wave function for visualization.

### Adjustable Parameters

- `L`: Domain size of the simulation grid.
- `Dy`: Grid spacing.
- `Dt`: Time step size.
- `v0`: Potential strength for the barrier.
- `x0, y0`: Initial position of the Gaussian wave packet.
- `w, s, a`: Geometry of the double-slit barrier.

## Running the Simulation

1. **Run the Script**: Execute the Python script to start the simulation.
2. **View the Animation**: The wave function's modulus will be animated, showing the propagation and interference.

## Example Output

The simulation visualizes the wave function evolving through the double-slit barrier, resulting in a pattern reminiscent of quantum interference. The resulting animation demonstrates how quantum mechanics predicts wave-like behavior in particle dynamics.

## Potential Applications

- Educational tool for visualizing quantum mechanical phenomena.
- Demonstration of sparse matrix methods in computational physics.
- Extension for simulating different potential configurations or multi-slit experiments.

## Extending the Project

- **Add More Complex Potentials**: Include time-dependent or irregular potentials.
- **Enhance Visualization**: Overlay plots showing the real and imaginary parts of the wave function.
- **Optimize Performance**: Parallelize the sparse matrix computations for larger grids.

## Acknowledgments

This project showcases the beauty of computational physics and its ability to simulate intricate quantum phenomena. Inspired by foundational experiments in quantum mechanics, it aims to bring abstract concepts to life through visualization.
