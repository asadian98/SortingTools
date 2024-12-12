# SortingTools

Welcome to the **SortingTools** repository! This project contains tools for **data extraction**, **spike sorting**, and **basic visualization** for neuroscience research. The code is designed to work seamlessly with data from **MonkeyLogic** and the **Ripple recording system**. 

## Overview

This repository provides:
- Scripts for extracting data from MonkeyLogic and Ripple.
- Spike sorting tools using popular external libraries.
- Basic analysis and visualization functions with an easy-to-use GUI.

The core pipeline is defined in **`KiloPipeline_main.m`**, with detailed comments about the outputs and file organization.

---

## Features

### 1. Spike Sorting
The pipeline integrates with several widely used spike sorting and quality assessment tools.

### 2. Data Analysis
Analysis and visualization functions are provided in the **`~\AnalysisUtils`** folder, including a demo script (`AnalysisUtils_demo.m`) with examples for each function.

### 3. Customizable
Easily extend the pipeline to include your own analysis tools and adapt it for different data sources.

---

## Dependencies

SortingTools leverages several external repositories for spike sorting and analysis. These are included in the repository for convenience:

| Tool              | Description                                    | Link                                                                                         |
|-------------------|------------------------------------------------|----------------------------------------------------------------------------------------------|
| **Kilosort**      | Fast and efficient spike sorting (supports Kilosort v2.0 and v4.0)               | [Kilosort v2.0](https://github.com/MouseLand/Kilosort/releases/tag/v2.0) [Kilosort version 4.0.22.dev5+g40059fc] (https://github.com/MouseLand/Kilosort)|
| **npy-matlab**    | Python-compatible .npy file reader for MATLAB  | [npy-matlab](https://github.com/kwikteam/npy-matlab)                                         |
| **Phy**           | GUI for manual curation of spike sorting       | [Phy](https://github.com/kwikteam/phy)                                                      |
| **spikes**        | Spike processing tools from Cortex Lab         | [spikes](https://github.com/cortex-lab/spikes)                                               |
| **sortingQuality**| Spike sorting quality metrics                 | [sortingQuality](https://github.com/cortex-lab/sortingQuality)                               |
| **kilo2Tools**    | Kilosort helper functions                      | [kilo2Tools](https://github.com/ElKatz/kilo2Tools)                                           |
| **CSDplotter**    | Current source density (CSD) visualization     | [CSDplotter](https://github.com/espenhgn/CSDplotter)                                         |
| **FastICA**       | Independent component analysis                 | [FastICA](https://github.com/aludnam/MATLAB/tree/master/FastICA_25)                          |
| **Bombcell**      | Spike sorting and analysis utilities           | [Bombcell](https://github.com/Julie-Fabre/bombcell/tree/main)                                |

---

## Getting Started

### Setup
1. Clone this repository to your local machine:
   ```bash
   git clone <repository-url>

2. Open MATLAB and navigate to the repository folder.
3. Run the demo script to explore the basic functionality:
    ```bash
    AnalysisUtils_demo.m
4. Start the pipeline with:
    ```bash
    KiloPipeline_main.m
### File Organization
Output files and their organization are explained in the comments within KiloPipeline_main.m.

### License
This project is open-source and distributed under the MIT License. Please review the license file for details.
