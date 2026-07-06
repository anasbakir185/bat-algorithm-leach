# Enhancement of LEACH Protocol in WSNs using the Bat Algorithm (BA-LEACH)

Graduation project, Communication Engineering Department, College of Electronic Technology, Tripoli, Libya.
Spring 2026.

**Student:** Anas Fauzi Bakeer
**Supervisor:** Dr. Ebtisam Mohamed Elgdiri

## Abstract

Wireless Sensor Networks consist of battery-powered sensors that may be deployed in large numbers to sense, measure, or gather information from the environment. A major challenge in deploying WSNs is the energy limitation of nodes. To address this challenge, energy efficient routing protocols such as LEACH are implemented to reduce the energy consumption of nodes during transmission. LEACH can achieve a significant increase in network lifetime, but it may suffer from early node failures due to its random cluster head selection process, which can elect nodes with low residual energy or unfavourable locations.

This project proposes BA-LEACH, an enhanced protocol that integrates the Bat Algorithm to select more suitable cluster heads each round based on residual energy and location. This approach mitigates uneven energy consumption and keeps more nodes functional for a longer time. BA-LEACH also introduces a multi-hop communication variant that reduces long-distance transmissions by allowing distant cluster heads to forward data through intermediate nodes, lowering overall energy consumption.

Simulations conducted in MATLAB show that single-hop BA-LEACH improves the stability period by 62%, while multi-hop BA-LEACH achieves a 94% stability improvement over standard LEACH.

## Key Results

| Protocol Variant   | Stability Period Improvement over LEACH |
|---------------------|------------------------------------------|
| BA-LEACH (Single-hop) | +62% |
| BA-LEACH (Multi-hop)  | +94% |

## Repository Structure

```
.
├── Code/
│   ├── LEACH.m           # Original LEACH protocol implementation
│   ├── BA-LEACH-S.m      # Bat Algorithm integrated into LEACH (single-hop)
│   ├── BA-LEACH-M.m      # Bat Algorithm integrated into LEACH (multi-hop)
│   ├── dist.m            # Distance calculation between two nodes
│   └── txEnergy.m        # Transmission energy calculation for one packet
└── README.md
```

## Prerequisites

- MATLAB (developed and tested with MATLAB R20XX, update to your version)

## How to Run

1. Clone or download this repository.
2. Open the `Code/` folder in MATLAB.
3. Run the scripts in order to compare protocol performance:
   1. `LEACH.m`
   2. `BA-LEACH-S.m`
   3. `BA-LEACH-M.m`
4. Results are displayed as plots and printed to the console.

## Notes

- All MATLAB code was written by hand, without the use of AI tools or code generation models.
- The same network topology and identical LEACH decision logic are used across all scenarios, ensuring a fair and consistent comparison of results.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Contact

- Email: anasbakir184@gmail.com
