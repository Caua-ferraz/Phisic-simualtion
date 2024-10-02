import tkinter as tk
from tkinter import messagebox
import math
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np

# Material properties (coefficients of restitution)
MATERIALS = {
    "Metal": 0.6,         # Less elastic
    "Rubber": 0.85,       # Highly elastic
    "Soft Plastic": 0.4   # Less elastic
}

# Direction vectors for launching the ball
DIRECTIONS = {
    "Up": (0, -1),
    "Down": (0, 1),
    "Left": (-1, 0),
    "Right": (1, 0)
}

# Gravitational accelerations for different celestial bodies (m/s²)
GRAVITIES = {
    "Earth": 9.81,
    "Moon": 1.62,
    "Mars": 3.71
}

# Air densities for different celestial bodies (kg/m³)
AIR_DENSITIES = {
    "Earth": 1.225,  # kg/m³ at sea level
    "Moon": 0.0,     # No atmosphere
    "Mars": 0.020    # kg/m³
}

# Colors associated with each celestial body for visualization
COLORS = {
    "Earth": "green",
    "Mars": "red",
    "Moon": "blue"
}

class Vector2D:
    """A simple 2D vector class for basic vector operations."""
    def __init__(self, x: float, y: float):
        self.x = x
        self.y = y

    def __add__(self, other: 'Vector2D') -> 'Vector2D':
        return Vector2D(self.x + other.x, self.y + other.y)

    def __sub__(self, other: 'Vector2D') -> 'Vector2D':
        return Vector2D(self.x - other.x, self.y - other.y)

    def __mul__(self, scalar: float) -> 'Vector2D':
        return Vector2D(self.x * scalar, self.y * scalar)

    def magnitude(self) -> float:
        return math.hypot(self.x, self.y)

    def normalized(self) -> 'Vector2D':
        mag = self.magnitude()
        if mag == 0:
            return Vector2D(0, 0)
        return Vector2D(self.x / mag, self.y / mag)

    def __str__(self) -> str:
        return f"Vector2D(x={self.x}, y={self.y})"

class Ball:
    """Represents the ball in the simulation with its physical properties and behaviors."""
    MIN_VELOCITY_THRESHOLD = 0.1  # m/s

    def __init__(self, gravity_name: str, initial_force: float, direction_vector: Vector2D,
                 mass: float, ball_radius_cm: float, drag_coefficient: float,
                 material_elasticity: float, scale: float, canvas: tk.Canvas,
                 room_size_m: float, air_density: float, canvas_height: float):
        self.gravity_name = gravity_name
        self.gravity = Vector2D(0, GRAVITIES[gravity_name])
        self.mass = mass
        self.ball_radius_m = ball_radius_cm / 100  # Convert cm to meters
        self.drag_coefficient = drag_coefficient
        self.air_density = air_density
        self.scale = scale
        self.canvas = canvas
        self.color = COLORS[gravity_name]
        self.material_elasticity = material_elasticity
        self.friction = 0.99  # Friction factor after each bounce
        self.room_size_m = room_size_m
        self.angular_velocity = Vector2D(0, 0)  # Simplified angular velocity
        self.moment_of_inertia = (2/5) * self.mass * (self.ball_radius_m ** 2)

        # Initialize position at the center of the room
        self.position = Vector2D(room_size_m / 2, room_size_m / 2)

        # Calculate initial velocity based on force and direction
        acceleration = direction_vector * (initial_force / self.mass)
        self.velocity = acceleration

        # Create the ball on the canvas
        self.ball_radius_px = self.ball_radius_m * self.scale
        self.canvas_height = canvas_height
        self.canvas_ball = self.canvas.create_oval(
            self.position.x * self.scale - self.ball_radius_px,
            self.invert_y(self.position.y * self.scale - self.ball_radius_px),
            self.position.x * self.scale + self.ball_radius_px,
            self.invert_y(self.position.y * self.scale + self.ball_radius_px),
            fill=self.color
        )

        # Data for plotting
        self.time_steps = []
        self.positions_y = []
        self.velocities_y = []
        self.speeds = []

        self.total_distance = 0.0
        self.bounce_count = 0
        self.is_moving = True

    def invert_y(self, y_coord: float) -> float:
        """Invert the y-coordinate for correct rendering on the Tkinter canvas."""
        return self.canvas_height - y_coord

    def compute_air_resistance(self) -> Vector2D:
        """Compute the air resistance acceleration based on current velocity."""
        area = math.pi * (self.ball_radius_m ** 2)
        speed = self.velocity.magnitude()
        if speed == 0:
            return Vector2D(0, 0)
        drag_force_magnitude = 0.5 * self.drag_coefficient * self.air_density * area * (speed ** 2)
        drag_force = self.velocity.normalized() * -drag_force_magnitude
        return drag_force * (1 / self.mass)

    def compute_magnus_effect(self) -> Vector2D:
        """Compute the Magnus effect (lift force due to spin)."""
        lift_coefficient = 0.1  # Tunable parameter
        angular_speed = self.angular_velocity.magnitude()
        speed = self.velocity.magnitude()
        if angular_speed == 0 or speed == 0:
            return Vector2D(0, 0)
        
        spin_parameter = (angular_speed * self.ball_radius_m) / speed
        lift_force_magnitude = 0.5 * self.air_density * (speed ** 2) * (math.pi * self.ball_radius_m ** 2) * lift_coefficient * spin_parameter
        lift_direction = Vector2D(-self.velocity.y, self.velocity.x).normalized()
        return lift_direction * (lift_force_magnitude / self.mass)

    def compute_wind_effect(self, wind_velocity: Vector2D) -> Vector2D:
        """Compute the acceleration due to wind."""
        relative_velocity = self.velocity - wind_velocity
        speed = relative_velocity.magnitude()
        if speed == 0:
            return Vector2D(0, 0)
        drag_force_magnitude = 0.5 * self.drag_coefficient * self.air_density * (math.pi * self.ball_radius_m ** 2) * (speed ** 2)
        drag_force = relative_velocity.normalized() * -drag_force_magnitude
        return drag_force * (1 / self.mass)

    def update(self, time_step: float, wind_velocity: Vector2D):
        """Update the ball's position and velocity based on forces acting on it."""
        if not self.is_moving:
            return

        # Calculate all accelerations
        gravity_acc = self.gravity
        air_resistance_acc = self.compute_air_resistance()
        magnus_acc = self.compute_magnus_effect()
        wind_acc = self.compute_wind_effect(wind_velocity)
        net_acc = gravity_acc + air_resistance_acc + magnus_acc + wind_acc

        # Update velocity and position
        self.velocity = self.velocity + (net_acc * time_step)
        displacement = self.velocity * time_step
        new_position = self.position + displacement

        # Update total distance traveled
        self.total_distance += displacement.magnitude()

        # Update angular velocity with decay
        self.angular_velocity = self.angular_velocity * 0.99  # Simple decay model

        # Update position
        self.position = new_position

        # Move the ball on the canvas with inverted y-coordinate
        self.canvas.coords(
            self.canvas_ball,
            self.position.x * self.scale - self.ball_radius_px,
            self.invert_y(self.position.y * self.scale - self.ball_radius_px),
            self.position.x * self.scale + self.ball_radius_px,
            self.invert_y(self.position.y * self.scale + self.ball_radius_px)
        )

        # Handle collisions with walls
        self.handle_collisions()

    def handle_collisions(self):
        """Detect and respond to collisions with the room walls."""
        room = self.room_size_m
        radius = self.ball_radius_m

        # Bottom wall
        if self.position.y + radius >= room:
            if self.velocity.y > 0:
                self.velocity.y = -abs(self.velocity.y) * self.material_elasticity
                self.velocity.x *= self.friction
                self.position.y = room - radius
                self.bounce_count += 1

        # Top wall
        if self.position.y - radius <= 0:
            if self.velocity.y < 0:
                self.velocity.y = abs(self.velocity.y) * self.material_elasticity
                self.position.y = radius

        # Right wall
        if self.position.x + radius >= room:
            if self.velocity.x > 0:
                self.velocity.x = -abs(self.velocity.x) * self.material_elasticity
                self.velocity.y *= self.friction
                self.position.x = room - radius
                self.bounce_count += 1

        # Left wall
        if self.position.x - radius <= 0:
            if self.velocity.x < 0:
                self.velocity.x = abs(self.velocity.x) * self.material_elasticity
                self.velocity.y *= self.friction
                self.position.x = radius
                self.bounce_count += 1

        # Stopping condition
        if (abs(self.velocity.x) < self.MIN_VELOCITY_THRESHOLD and
            abs(self.velocity.y) < self.MIN_VELOCITY_THRESHOLD and
            self.position.y + radius >= room - 0.01):
            self.is_moving = False
            self.position.y = room - radius
            # Update the ball's position on the canvas one last time
            self.canvas.coords(
                self.canvas_ball,
                self.position.x * self.scale - self.ball_radius_px,
                self.invert_y(self.position.y * self.scale - self.ball_radius_px),
                self.position.x * self.scale + self.ball_radius_px,
                self.invert_y(self.position.y * self.scale + self.ball_radius_px)
            )
        else:
            pass  # Continue moving

    def record_data(self, current_time: float):
        """Record data for plotting."""
        self.time_steps.append(current_time)
        self.positions_y.append(self.position.y)
        self.velocities_y.append(self.velocity.y)
        self.speeds.append(self.velocity.magnitude() * 3.6)  # Convert m/s to km/h

class PhysicsSimulation:
    """Main class to handle the physics simulation application."""
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("Ball Launch Physics Simulation")
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)

        # Simulation parameters
        self.canvas_width = 600
        self.canvas_height = 400
        self.time_step = 0.02  # 50 frames per second

        # Initialize components
        self.balls = []
        self.time_elapsed = 0.0
        self.is_moving = False

        # Create GUI components
        self.create_controls()
        self.create_canvas()
        self.create_plots()

    def create_controls(self):
        """Create the control panel for simulation inputs."""
        control_frame = tk.Frame(self.root)
        control_frame.pack(side=tk.LEFT, padx=10, pady=10, fill=tk.Y)

        # Material selection
        tk.Label(control_frame, text="Material:").grid(row=0, column=0, padx=5, pady=5, sticky='w')
        self.material_var = tk.StringVar(value="Rubber")
        material_menu = tk.OptionMenu(control_frame, self.material_var, *MATERIALS.keys())
        material_menu.grid(row=0, column=1, padx=5, pady=5)

        # Launch force
        tk.Label(control_frame, text="Launch Force (N):").grid(row=1, column=0, padx=5, pady=5, sticky='w')
        self.force_var = tk.DoubleVar(value=10.0)
        tk.Entry(control_frame, textvariable=self.force_var).grid(row=1, column=1, padx=5, pady=5)

        # Launch direction
        tk.Label(control_frame, text="Launch Direction:").grid(row=2, column=0, padx=5, pady=5, sticky='w')
        self.direction_var = tk.StringVar(value="Up")
        direction_menu = tk.OptionMenu(control_frame, self.direction_var, *DIRECTIONS.keys())
        direction_menu.grid(row=2, column=1, padx=5, pady=5)

        # Room size
        tk.Label(control_frame, text="Room Size (m):").grid(row=3, column=0, padx=5, pady=5, sticky='w')
        self.room_size_var = tk.DoubleVar(value=10.0)
        tk.Entry(control_frame, textvariable=self.room_size_var).grid(row=3, column=1, padx=5, pady=5)

        # Ball size
        tk.Label(control_frame, text="Ball Size (cm):").grid(row=4, column=0, padx=5, pady=5, sticky='w')
        self.ball_size_var = tk.DoubleVar(value=5.0)
        tk.Entry(control_frame, textvariable=self.ball_size_var).grid(row=4, column=1, padx=5, pady=5)

        # Mass
        tk.Label(control_frame, text="Mass (kg):").grid(row=5, column=0, padx=5, pady=5, sticky='w')
        self.mass_var = tk.DoubleVar(value=1.0)
        tk.Entry(control_frame, textvariable=self.mass_var).grid(row=5, column=1, padx=5, pady=5)

        # Drag coefficient
        tk.Label(control_frame, text="Drag Coefficient:").grid(row=6, column=0, padx=5, pady=5, sticky='w')
        self.drag_coefficient_var = tk.DoubleVar(value=0.47)  # Default for a sphere
        tk.Entry(control_frame, textvariable=self.drag_coefficient_var).grid(row=6, column=1, padx=5, pady=5)

        # Wind speed input
        tk.Label(control_frame, text="Wind Speed (m/s):").grid(row=7, column=0, padx=5, pady=5, sticky='w')
        self.wind_speed_var = tk.DoubleVar(value=0.0)
        tk.Entry(control_frame, textvariable=self.wind_speed_var).grid(row=7, column=1, padx=5, pady=5)

        # Wind direction selection
        tk.Label(control_frame, text="Wind Direction:").grid(row=8, column=0, padx=5, pady=5, sticky='w')
        self.wind_direction_var = tk.StringVar(value="Right")
        wind_direction_menu = tk.OptionMenu(control_frame, self.wind_direction_var, "Left", "Right")
        wind_direction_menu.grid(row=8, column=1, padx=5, pady=5)

        # Temperature input
        tk.Label(control_frame, text="Temperature (°C):").grid(row=9, column=0, padx=5, pady=5, sticky='w')
        self.temperature_var = tk.DoubleVar(value=20.0)
        tk.Entry(control_frame, textvariable=self.temperature_var).grid(row=9, column=1, padx=5, pady=5)

        # Humidity input
        tk.Label(control_frame, text="Relative Humidity (%):").grid(row=10, column=0, padx=5, pady=5, sticky='w')
        self.humidity_var = tk.DoubleVar(value=50.0)
        tk.Entry(control_frame, textvariable=self.humidity_var).grid(row=10, column=1, padx=5, pady=5)

        # Throw Ball Button
        tk.Button(control_frame, text="Throw Ball", command=self.start_simulation).grid(row=11, column=0, columnspan=2, pady=10)

    def create_canvas(self):
        """Create the simulation canvas where the ball moves."""
        self.canvas = tk.Canvas(self.root, width=self.canvas_width, height=self.canvas_height, bg='white')
        self.canvas.pack(side=tk.RIGHT, padx=10, pady=10)

    def create_plots(self):
        """Create real-time plots for the simulation data."""
        fig, (self.ax1, self.ax2, self.ax3) = plt.subplots(3, 1, figsize=(5, 8))
        plt.tight_layout(pad=3.0)

        self.canvas_fig = FigureCanvasTkAgg(fig, master=self.root)
        self.canvas_fig.draw()
        self.canvas_fig.get_tk_widget().pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        self.fig = fig  # Reference to the figure for updating

    def calculate_air_density(self, temperature: float, humidity: float) -> float:
        """Calculate air density based on temperature and humidity."""
        pressure = 101325  # Standard atmospheric pressure in Pa
        R_dry = 287.05     # Specific gas constant for dry air in J/(kg·K)
        temperature_k = temperature + 273.15  # Convert Celsius to Kelvin

        # Simplified saturation vapor pressure (Tetens formula)
        e_s = 610.78 * math.exp((17.27 * temperature) / (temperature + 237.3))
        e = (humidity / 100) * e_s

        # Air density calculation
        density = (pressure - 0.378 * e) / (R_dry * temperature_k)
        return density

    def start_simulation(self):
        """Initialize and start the simulation."""
        if self.is_moving:
            messagebox.showinfo("Simulation Running", "The simulation is already running.")
            return

        try:
            # Gather input parameters
            material = self.material_var.get()
            initial_force = self.force_var.get()
            direction = self.direction_var.get()
            room_size = self.room_size_var.get()
            ball_size = self.ball_size_var.get()
            mass = self.mass_var.get()
            drag_coefficient = self.drag_coefficient_var.get()
            wind_speed = self.wind_speed_var.get()
            wind_direction = self.wind_direction_var.get()
            temperature = self.temperature_var.get()
            humidity = self.humidity_var.get()

            # Validate inputs
            if mass <= 0 or ball_size <= 0 or room_size <= 0:
                messagebox.showerror("Invalid Input", "Mass, Ball Size, and Room Size must be positive numbers.")
                return

            # Calculate air density
            air_density = self.calculate_air_density(temperature, humidity)

            # Determine wind velocity vector
            wind_multiplier = 1 if wind_direction == "Right" else -1
            wind_velocity = Vector2D(wind_speed * wind_multiplier, 0)

            # Determine launch direction vector
            direction_vector = Vector2D(*DIRECTIONS[direction])

            # Calculate scale to fit the room within the canvas
            scale = min(self.canvas_width / room_size, self.canvas_height / room_size)

            # Clear existing balls if any
            self.clear_balls()

            # Create a ball for each celestial body
            for gravity_name in GRAVITIES.keys():
                ball = Ball(
                    gravity_name=gravity_name,
                    initial_force=initial_force,
                    direction_vector=direction_vector,
                    mass=mass,
                    ball_radius_cm=ball_size,
                    drag_coefficient=drag_coefficient,
                    material_elasticity=MATERIALS.get(material, 0.5),
                    scale=scale,
                    canvas=self.canvas,
                    room_size_m=room_size,
                    air_density=AIR_DENSITIES.get(gravity_name, air_density),
                    canvas_height=self.canvas_height
                )
                self.balls.append(ball)

            self.is_moving = True
            self.time_elapsed = 0.0
            self.animate_balls(wind_velocity)

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {e}")

    def clear_balls(self):
        """Remove all balls from the canvas and reset the list."""
        for ball in self.balls:
            self.canvas.delete(ball.canvas_ball)
        self.balls.clear()

    def animate_balls(self, wind_velocity: Vector2D):
        """Animate the balls by updating their positions and plotting data."""
        if not self.is_moving:
            return

        all_stopped = True

        for ball in self.balls:
            if ball.is_moving:
                ball.update(self.time_step, wind_velocity)
                ball.record_data(self.time_elapsed)
                all_stopped = False

        # Update plots with the latest data
        self.update_plots()

        self.time_elapsed += self.time_step

        if all_stopped:
            self.is_moving = False
            self.display_results()
        else:
            self.root.after(int(self.time_step * 1000), lambda: self.animate_balls(wind_velocity))

    def update_plots(self):
        """Update the real-time plots with data from all balls."""
        # Clear previous plots
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()

        # Configure each subplot
        self.ax1.set_title('Vertical Position over Time')
        self.ax1.set_xlabel('Time (s)')
        self.ax1.set_ylabel('Position (m)')
        self.ax1.invert_yaxis()  # Invert to match the simulation's y-axis

        self.ax2.set_title('Vertical Velocity over Time')
        self.ax2.set_xlabel('Time (s)')
        self.ax2.set_ylabel('Velocity (m/s)')

        self.ax3.set_title('Speed over Time')
        self.ax3.set_xlabel('Time (s)')
        self.ax3.set_ylabel('Speed (km/h)')

        # Plot data for each ball
        for ball in self.balls:
            self.ax1.plot(ball.time_steps, ball.positions_y, label=ball.gravity_name, color=ball.color)
            self.ax2.plot(ball.time_steps, ball.velocities_y, label=ball.gravity_name, color=ball.color)
            self.ax3.plot(ball.time_steps, ball.speeds, label=ball.gravity_name, color=ball.color)

        # Add legends
        self.ax1.legend()
        self.ax2.legend()
        self.ax3.legend()

        # Redraw the canvas with updated plots
        self.fig.canvas.draw()

    def display_results(self):
        """Display the simulation results once all balls have stopped."""
        y_position = self.canvas_height / 2 - 20
        for ball in self.balls:
            result_text = (f"{ball.gravity_name}: Stopped after {ball.bounce_count} bounces, "
                           f"Distance: {ball.total_distance:.2f} m")
            self.canvas.create_text(
                self.canvas_width / 2,
                y_position,
                text=result_text,
                fill=ball.color,
                font=("Arial", 12)
            )
            y_position += 20

    def on_closing(self):
        """Handle the window closing event."""
        if self.is_moving:
            if messagebox.askokcancel("Quit", "Simulation is running. Do you want to quit?"):
                self.is_moving = False
                self.root.destroy()
        else:
            self.root.destroy()

def main():
    """Initialize and run the physics simulation application."""
    root = tk.Tk()
    app = PhysicsSimulation(root)
    root.mainloop()

if __name__ == "__main__":
    main()