# Ball Launch Physics Simulation with Walls

## Table of Contents

- [Introduction](#introduction)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Running the Simulation](#running-the-simulation)
- [User Interface](#user-interface)
  - [Simulation Parameters](#simulation-parameters)
  - [Controls](#controls)
- [Technical Description](#technical-description)
  - [Units and Scale](#units-and-scale)
  - [Physics Calculations](#physics-calculations)
    - [Gravity](#gravity)
    - [Air Resistance](#air-resistance)
    - [Magnus Effect](#magnus-effect)
    - [Wind Effect](#wind-effect)
    - [Collisions](#collisions)
  - [Animation and Updates](#animation-and-updates)
- [Graphs and Visualization](#graphs-and-visualization)
- [Customization and Extension](#customization-and-extension)
- [Troubleshooting](#troubleshooting)
- [References](#references)
- [License](#license)

## Introduction

Imagine throwing a ball in a room and watching it bounce around. Now, imagine doing that on Earth, Mars, and the Moon all at once! That's what our physics simulation lets you do. It's like having a virtual playground where you can experiment with different forces and see how they affect the ball's movement.

This simulation brings to life important concepts from classical physics, including:

- Two-dimensional (2D) motion
- Gravity's pull
- Air resistance (drag)
- The Magnus effect (spin-induced lift)
- Wind effects
- Bouncy and not-so-bouncy collisions with walls
- How mass and applied force change things

Our goal is to create an educational tool that makes understanding physics principles fun and interactive. It's like a physics lab in your computer!

## System Requirements

To run this simulation, you'll need:

- **Python 3.x**: The simulation is written in Python, so you'll need version 3 or higher.
- **Python Libraries**:
  - `tkinter`: For creating the graphical interface (usually comes with Python).
  - `matplotlib`: For creating those cool graphs during the simulation.
  - `numpy`: For some fancy number crunching.

Make sure you have these libraries installed and ready to go!

## Installation

Getting started is easy:

1. **Get the code**: Clone or download the source code from the repository.

2. **Install the dependencies**:
   Open your command prompt or terminal and run:

   ```bash
   pip install -r requirements.txt
   ```

   This will install all the necessary libraries listed in the `requirements.txt` file.

   *Note*: If you're using Linux, you might need to install `tkinter` separately.

3. **Run the script**:
   In the same command prompt or terminal, type:

   ```bash
   python main.py
   ```

   And watch the magic happen!

## Running the Simulation

Once you run the script, a window will pop up with all the controls on the left and the simulation area on the right. It's time to play physicist!

## User Interface

The interface is split into two main parts:

1. **Controls**: On the left, where you can tweak all the simulation settings.
2. **Simulation Area**: On the right, where you'll see the ball(s) bouncing around.

### Simulation Parameters

Before you launch the ball, you can adjust these settings:

- **Material**: Choose what the ball is made of (affects bounciness).
  - Options: *Rubber*, *Metal*, *Soft Plastic*
- **Launch Force (N)**: How hard you're throwing the ball (in Newtons).
- **Launch Direction**: Which way you're throwing.
  - Options: *Up*, *Down*, *Left*, *Right*
- **Room Size (m)**: How big is your virtual room (in meters).
- **Ball Size (cm)**: How big is your ball (diameter in centimeters).
- **Mass (kg)**: How heavy is your ball (in kilograms).
- **Drag Coefficient**: How much air resistance affects the ball.
- **Wind Speed (m/s)**: How fast is the wind blowing.
- **Wind Direction**: Which way is the wind blowing.
  - Options: *Left*, *Right*
- **Temperature (°C)**: How hot or cold is it (affects air density).
- **Relative Humidity (%)**: How muggy is it (also affects air density).

### Controls

- **"Throw Ball" Button**: Click this to start the show!

## Technical Description

### Units and Scale

We're keeping it scientific:

- **Units**: We use the International System of Units (SI) for all calculations.
  - Distance: meters (m)
  - Mass: kilograms (kg)
  - Time: seconds (s)
  - Force: Newtons (N)
  - Speed: meters per second (m/s)
  - Acceleration: meters per second squared (m/s²)
- **Visual Scale**: The simulation adjusts automatically so you can always see the whole room and ball, no matter their actual size.

### Physics Calculations

#### Gravity

- We use real gravity values for different celestial bodies:
  - Earth: 9.81 m/s²
  - Moon: 1.62 m/s²
  - Mars: 3.71 m/s²
  
  Gravity always pulls down (positive y-direction).

#### Air Resistance

- We calculate drag force using this equation:
  
  F_d = 0.5 * C_d * ρ * A * v²
  
  Where:
  - C_d is the drag coefficient
  - ρ is air density
  - A is the ball's cross-sectional area
  - v is the ball's speed
  
- The resulting deceleration opposes the ball's motion.

#### Magnus Effect

- This is the lift force due to the ball's spin:
  
  F_m = 0.5 * C_l * ρ * A * v² * ω
  
  Where:
  - C_l is the lift coefficient
  - ω is the ball's angular velocity
  
  This force can make the ball curve in flight.

#### Wind Effect

- We calculate how wind pushes the ball:
  
  F_w = 0.5 * C_d * ρ * A * (v - v_wind)²
  
  This can significantly affect the ball's trajectory.

#### Collisions

- **Bounciness**: Depends on the material you chose.
- **Wall Collisions**: The ball bounces off walls, losing some energy each time.
- **Friction**: We apply a friction factor to slow the ball down over time.

### Animation and Updates

- **Time Steps**: The simulation updates every 0.02 seconds (50 frames per second).
- **Position Updates**: We use Euler's method to calculate the ball's new position each step.
- **Visual Updates**: The ball's position on screen is updated each step.

## Graphs and Visualization

We provide real-time graphs to help you understand the ball's motion:

- **Position vs. Time**: Shows how high the ball is over time.
- **Velocity vs. Time**: Shows how fast the ball is moving vertically.
- **Speed vs. Time**: Shows the ball's overall speed in km/h.

These graphs update in real-time as the simulation runs.

## Customization and Extension

Want to take it further? Here are some ideas:

- **New Materials**: Add new ball materials with different bounciness levels.
- **Different Environments**: Adjust gravity, air resistance, etc., to simulate other planets or conditions.
- **More Features**:
  - Add pause and resume buttons.
  - Simulate multiple interacting balls.
  - Show graphs of kinetic and potential energy.
  - Add a feature to save simulation data for later analysis.

## Troubleshooting

Having issues? Check these common problems:

- **Ball not moving**: Make sure all your input values are valid (no zero mass, negative sizes, etc.).
- **Ball disappearing**: Check if your room size settings are reasonable.
- **Simulation won't start**: Ensure all required libraries are installed correctly.
- **Graphs not updating**: Check if your matplotlib backend is configured correctly.

## References

Want to learn more? Check out these resources:

- [Tkinter Documentation](https://docs.python.org/3/library/tkinter.html)
- [Matplotlib: Python plotting](https://matplotlib.org/)
- [NumPy Documentation](https://numpy.org/doc/)
- Physics Concepts:
  - Classical Mechanics: Study of object motion under forces.
  - Fluid Dynamics: Study of fluids in motion, including air resistance.
  - Magnus Effect: Lift force on spinning objects in motion.

## License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details. Feel free to use, modify, and share!