 # Real-Time Structural Health Monitoring of an Offshore Wind Turbine Using a Digital Twin

 ## Background and Motivation

 Offshore wind energy has rapidly emerged as a cornerstone of global decarbonization efforts. As offshore wind farms scale up in size, complexity, and remoteness, ensuring the reliability and structural integrity of wind turbines becomes increasingly critical. Real-time monitoring of turbine health, through Structural Health Monitoring (SHM) systems, plays a vital role in minimizing maintenance costs, predicting failures, and optimizing operational availability.

This personal project grew out of my master's thesis, where I developed a digital twin of an offshore wind turbine (OWT). Now, I have extended it into a full end-to-end SHM simulation platform — aiming to emulate an industry-grade monitoring system using open-source tools and Python. The goal: bridge the gap between simulation and real-time condition monitoring.

## Objective

The core objective of the project is to:

- Simulate a realistic offshore wind turbine SHM system using a digital twin
- Stream simulated sensor data (acceleration, velocity, displacement) in real time
- Store and manage the data using a structured database
- Visualize and monitor this data via a dashboard resembling a SCADA interface
- Perform data analysis to detect abnormal trends or potential structural issues

Ultimately, this serves as a **proof-of-concept for condition-based monitoring and digital twin integration in offshore wind** which is a valuable stepping stone toward smart maintenance strategies in the renewable energy industry.

## Project Architecture Overview

The project is structured around five key components:

1. Digital Twin Modeling (Simulation Layer)
2. Real-Time Data Emulation (Streaming Layer)
3. Database Storage (Persistence Layer)
4. Dashboarding (Monitoring Layer)
5. Data Analysis (Insight Layer)

## 1. Digital Twin Modelling

This project began with a digital twin of an offshore wind turbine tower that I had developed for my thesis. The twin simulates the dynamic response of the turbine tower under initial excitation and operating conditions, outputting:

- Acceleration
- Velocity
- Displacement

across a predefined time period, based on physical equations of motion. The simulation is implemented in Python using numerical solvers, and the output mimics what real accelerometers or vibration sensors would record at various locations on the turbine structure.

## 2. Real-Time Data Streaming

Rather than processing all the simulation data at once, the project is designed to simulate live sensor behavior. This is achieved by feeding the digital twin output into a custom Python script that:

- Loops through the simulated time-series data point by point
- Inserts each data point into a connected database at intervals (e.g., 10 Hz)
- Mimics the behavior of a live offshore sensor system

To make the system realistic and testable, I added logic for:

- Timestamping
- Simulating anomalies (e.g., sudden vibration spikes)
- Tagging data streams based on sensor ID or location

This streaming simulation forms the backbone of the “sensor-to-database” pathway.

## 3. Data storage

The next challenge is selecting an appropriate database for storing the streamed data. Two approaches were evaluated:

- **MySQL:** A traditional relational database that I was already comfortable with. It supports structured schema (timestamp, acceleration, velocity, displacement) and allows basic querying for historical data.
- **InfluxDB:** A purpose-built time-series database optimized for sensor data, with built-in features like retention policies, downsampling, and native Grafana integration.

Currently, I'm planning to use MySQL for its simplicity, but in future as the project gets bigger, I have plans of migrating the database to InfluxDB as its built for storing time-series data.

## 4. Real-Time Dashboarding

To monitor the health of the turbine visually, I set up Grafana, a powerful open-source dashboarding tool that resembles the functionality of an industrial **SCADA (Supervisory Control and Data Acquisition)** system.

Grafana was configured to:

- Pull real-time data from InfluxDB
- Display acceleration, velocity, and displacement curves on live plots
- Trigger alerts when values exceed specified thresholds
- Allow exploration of historical data trends

The dashboard replicates a control room feel, giving a clear view of structural dynamics and anomalies in real time. This integration showcases how simulated sensor data can be effectively turned into meaningful operational insight.

## 5. Data Analysis for Health Insights

The final and most exciting component is data analysis. Once the data is stored, I perform post-processing using Python tools like:

- NumPy / SciPy for signal analysis (FFT, RMS, zero-crossing rates)
- Pandas for trend evaluation
- Matplotlib / Plotly for advanced visualization

Some of the key analyses I experimented with include:

- Identifying frequency content from vibration data (for detecting structural resonance)
- Detecting sudden peaks indicating potential damage
- Tracking long-term drift in displacement or velocity, which might suggest loosening bolts, corrosion, or fatigue damage

This analysis layer helps simulate a condition-based maintenance strategy which is the core goal of SHM.

## Industrial Relavence

What makes this project unique is how it bridges theory with industry reality:

- It replicates key components of actual SHM pipelines used in offshore wind turbines
- Uses tools and workflows similar to those adopted in digital O&M platforms
- Aligns closely with trends like digital twins, predictive maintenance, and data-driven wind farm optimization

The architecture mirrors how modern offshore operators are integrating SCADA data, monitoring dashboards, and AI-based diagnostics to improve turbine reliability while reducing O&M costs.