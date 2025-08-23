from simulation import TowerSimulation
from force_calculator import ForceCalculator

def main():
    simulation = TowerSimulation()
    force_calculator = ForceCalculator()
    print(f"x = {force_calculator.val}")
    print(f"f_wind = {force_calculator.force_wind()}")
    print(f"f_aero = {force_calculator.force_aero()}")
    print(f"f_wave = {simulation.force_wave_values}")
    # Run simulation
    simulation.run()

if __name__ == "__main__":
    main()