# Fluid Simulation - OpenMP

This project is developed as a part of the Semester 6 Heterogeneous Parallelism Course. It aims to parallelize certain parts of fluid dynamics simulations using OpenMP to improve computational efficiency.

## Problem Statement

Fluid dynamics simulations involve solving complex partial differential equations, such as the Navier-Stokes equations and the Gauss-Seidel model, which are computationally intensive. When performed serially, these simulations can take a long time due to high computational demand and dependencies.

## Solution

To address this issue, we have utilized OpenMP to parallelize specific segments of the simulation, allowing for concurrent execution on multiple threads and thus reducing the overall computational time.

## Demo

    ![here](demo/demo.gif).

## References

- [Real-time Fluid Dynamics for Games](https://damassets.autodesk.net/content/dam/autodesk/www/autodesk-reasearch/Publications/pdf/realtime-fluid-dynamics-for.pdf)
- [Fluid Simulation for Dummies](https://mikeash.com/pyblog/fluid-simulation-for-dummies.html)

## Contributors

- Gayathri Manoj
- Sai Manasa Nadimpalli