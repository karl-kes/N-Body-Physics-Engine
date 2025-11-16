#!/usr/bin/env python3
"""
Advanced N-Body Initial Conditions Generator
Generates various interesting gravitational scenarios for simulation
"""

import numpy as np
from typing import List, Tuple
import sys

class BodyGenerator:
    def __init__(self, seed=None):
        if seed is not None:
            np.random.seed(seed)
    
    def write_bodies(self, bodies: List[Tuple], filename: str, description: str):
        """Write bodies to CSV file with header"""
        with open(filename, 'w') as f:
            f.write(f"# {description}\n")
            for body in bodies:
                f.write(','.join(map(str, body)) + '\n')
        print(f"Generated {len(bodies)} bodies in '{filename}'")
    
    def planetary_system_with_rings(self, star_mass=2e30):
        """Generate a star with planets, including one with a ring system"""
        bodies = []
        G = 6.67430e-11
        
        # Central star
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Inner rocky planets
        for i in range(3):
            dist = (i + 1) * 5e10
            angle = np.random.uniform(0, 2*np.pi)
            mass = 10**np.random.uniform(24, 25)
            
            x = dist * np.cos(angle)
            y = dist * np.sin(angle)
            z = 0
            
            v = np.sqrt(G * star_mass / dist)
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        # Gas giant with rings (like Saturn)
        giant_dist = 1.5e11
        giant_mass = 5e27
        giant_x = giant_dist
        giant_y = 0
        giant_vx = 0
        giant_vy = np.sqrt(G * star_mass / giant_dist)
        
        bodies.append((giant_x, giant_y, 0, giant_vx, giant_vy, 0, giant_mass))
        
        # Create ring system (many small bodies)
        ring_inner = 7e7
        ring_outer = 1.5e8
        n_ring_particles = 100
        
        for i in range(n_ring_particles):
            r = np.random.uniform(ring_inner, ring_outer)
            theta = np.random.uniform(0, 2*np.pi)
            
            # Position relative to giant
            x = giant_x + r * np.cos(theta)
            y = giant_y + r * np.sin(theta)
            z = np.random.normal(0, 1e6)  # Very thin ring
            
            # Orbital velocity around giant + giant's velocity
            v_ring = np.sqrt(G * giant_mass / r)
            vx = giant_vx - v_ring * np.sin(theta)
            vy = giant_vy + v_ring * np.cos(theta)
            
            mass = 10**np.random.uniform(18, 20)  # Small particles
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        # Outer ice giants
        for i in range(2):
            dist = (3 + i) * 1e11
            angle = np.random.uniform(0, 2*np.pi)
            mass = 10**np.random.uniform(26, 27)
            
            x = dist * np.cos(angle)
            y = dist * np.sin(angle)
            z = np.random.uniform(-1e10, 1e10)
            
            v = np.sqrt(G * star_mass / dist)
            vx = -v * np.sin(angle)
            vy = v * np.cos(angle)
            
            bodies.append((x, y, z, vx, vy, 0, mass))
        
        return bodies
    
    def globular_cluster(self, n_bodies=1000, radius=1e13, core_radius=2e12):
        """Generate a globular cluster with core collapse dynamics"""
        bodies = []
        G = 6.67430e-11
        total_mass = n_bodies * 1e30  # Approximate
        
        for i in range(n_bodies):
            # King model-like distribution (denser core)
            if np.random.random() < 0.3:  # 30% in core
                r = core_radius * (np.random.random() ** (1/3))
            else:  # 70% in halo
                r = core_radius + (radius - core_radius) * (np.random.random() ** (1/2))
            
            # Random direction
            theta = np.arccos(2 * np.random.random() - 1)
            phi = 2 * np.pi * np.random.random()
            
            x = r * np.sin(theta) * np.cos(phi)
            y = r * np.sin(theta) * np.sin(phi)
            z = r * np.cos(theta)
            
            # Velocity dispersion decreases with radius
            sigma = np.sqrt(G * total_mass / (r + core_radius))
            vx = np.random.normal(0, sigma * 0.3)
            vy = np.random.normal(0, sigma * 0.3)
            vz = np.random.normal(0, sigma * 0.3)
            
            # Variable stellar masses (main sequence to giants)
            if i < 5:  # Few massive stars/black holes
                mass = 10**np.random.uniform(30, 31)
            elif i < 20:  # Giants
                mass = 10**np.random.uniform(29, 30)
            else:  # Main sequence
                mass = 10**np.random.uniform(28, 29)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies
    
    def protoplanetary_disk(self, n_particles=150):
        """Young star with protoplanetary disk that will form planets"""
        bodies = []
        G = 6.67430e-11
        star_mass = 2e30
        
        # Young star
        bodies.append((0, 0, 0, 0, 0, 0, star_mass))
        
        # Disk particles with Keplerian rotation
        inner_radius = 3e10
        outer_radius = 5e11
        
        for i in range(n_particles):
            # Power law distribution (more particles closer in)
            r = inner_radius * (outer_radius/inner_radius) ** np.random.random()
            theta = np.random.uniform(0, 2*np.pi)
            
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            # Thin disk with scale height proportional to radius
            z = np.random.normal(0, r * 0.01)
            
            # Keplerian velocity
            v_kep = np.sqrt(G * star_mass / r)
            # Add some radial drift
            v_r = np.random.normal(0, v_kep * 0.01)
            
            vx = -v_kep * np.sin(theta) + v_r * np.cos(theta)
            vy = v_kep * np.cos(theta) + v_r * np.sin(theta)
            vz = np.random.normal(0, v_kep * 0.001)
            
            # Mass distribution (dust to planetesimals)
            if i < 10:  # Few large planetesimals
                mass = 10**np.random.uniform(23, 25)
            else:  # Mostly dust and small particles
                mass = 10**np.random.uniform(18, 22)
            
            bodies.append((x, y, z, vx, vy, vz, mass))
        
        return bodies

def main():
    # Randomly pick a scenario
    gen = BodyGenerator()
    
    scenarios = {
        'rings': (gen.planetary_system_with_rings, "Planetary System with Ringed Giant"),
        'cluster': (gen.globular_cluster, "Globular Cluster"),
        'disk': (gen.protoplanetary_disk, "Protoplanetary Disk")
    }
    
    # Pick a random scenario
    scenario_name = np.random.choice(list(scenarios.keys()))
    func, desc = scenarios[scenario_name]
    
    print(f"Generating scenario: {desc}")
    print("=" * 50)
    
    # Generate the bodies
    bodies = func()
    
    # Write to bodies.csv
    gen.write_bodies(bodies, 'bodies.csv', desc)
    
    # Print some statistics
    masses = [b[6] for b in bodies]
    print(f"\nStatistics:")
    print(f"  Total bodies: {len(bodies)}")
    print(f"  Mass range: {min(masses):.2e} - {max(masses):.2e} kg")
    print(f"  Total mass: {sum(masses):.2e} kg")
    print(f"\nTip: Run again to generate a different random scenario!")

if __name__ == "__main__":
    main()