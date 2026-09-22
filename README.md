assemblyCount
A Python-based tool for computing genome assembly statistics


https://img.shields.io/badge/python-3.9%2B-blue.svg
https://img.shields.io/badge/License-MIT-green.svg
https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg
https://img.shields.io/badge/version-1.0-orange.svg
https://img.shields.io/badge/GUI-Tkinter-blue.svg


assemblyCount is a Python-based tool that evaluates the quality and completeness of de novo genome assemblies. It parses multi-FASTA files to compute essential assembly metrics (Total Size, N50, L50, N90, L90, and GC content) and generates a highly customizable, publication-ready cumulative scaffold bar chart.


The tool provides a modern, blue-themed graphical user interface with real-time logging, making it accessible to researchers without command-line expertise, while maintaining the rigor required for genomic analysis.


The pipeline combines:

Robust parsing of multi-FASTA files (large file tolerant)

Accurate calculation of Nx and Lx metrics (N50, L50, N90, L90)

GC content estimation

Publication-quality cumulative bar chart generation (300 DPI, PDF/SVG support)

Threaded GUI execution to prevent freezing on large genomes

Real-time GUI console with color-coded messages and progress bar

Detailed TXT report generation


Table of Contents
Overview

What  assemblyCount Does

New in Version 1.0

Pipeline Workflow

System Requirements

Installation Guide

Understanding Assembly Metrics

Input Files

Output Files

Using the Graphical Interface

Using the Command-Line Interface

Understanding the Results

Troubleshooting

FAQ

Citation

License

Version History



Overview
Genome assemblies produced by long-read or hybrid pipelines require rigorous quality assessment before downstream annotation. 
Key metrics such as N50 and L50 provide a measure of contiguity, while N90 and L90 offer insight into the completeness of the assembly. 
Additionally, visualizing the cumulative contribution of scaffolds helps researchers quickly identify the largest contigs and assess overall assembly fragmentation.



assemblyCount addresses this need in three steps:

Parses the input multi-FASTA and sorts scaffolds by size in descending order.

Computes Total Size, N50, L50, N90, L90, and GC content.

Generates a cumulative bar chart (Green = sum of previous scaffolds, Red = current scaffold) and exports it in high resolution.


The tool is organism-agnostic and requires only a standard FASTA file. No reference genome, no annotation, and no external alignment software.



What  assemblyCount Does
Core Functions
FASTA Parsing
Reads standard multi-FASTA files (.fasta, .fa, .fna).

Efficiently handles large genomes by reading line-by-line and storing sequence lengths.

Preserves scaffold identifiers and descriptions.


Assembly Metric Calculation
Sorts scaffolds by length (descending).

Computes Total Assembly Size.

Calculates N50 and L50 (the length and number of scaffolds at 50% of the assembly).

Calculates N90 and L90 (the length and number of scaffolds at 90% of the assembly).

Calculates overall GC content percentage.


Publication-Quality Visualization
Generates a cumulative bar chart using Matplotlib.

Green bars represent the cumulative sum of all previous scaffolds.

Red bars represent the current scaffold.

Limits plotting to the top 50 scaffolds by default to maintain readability.

Saves figures in PNG (300 DPI), PDF, or SVG (vector) formats.


Reporting
Generates a plain text report (.txt) containing all calculated statistics and the top 5 largest scaffolds.


Graphical Interface
Modern blue-themed Tkinter GUI.

Threaded execution: the interface remains responsive during large analyses.

Real-time color-coded log console.

Progress bar.

Scrollable form panel so every control stays reachable on small screens.


Command-Line Interface
(Optional) Can be scripted to run headless for pipeline integration.

New in Version 1.0
First public release of  assemblyCount.

Unified GUI and logic in a single portable script.

Robust calculation of N50, L50, N90, L90, and GC content.

Publication-quality cumulative bar chart generation.

Detailed TXT report.



Pipeline Workflow
text
INPUT FILES
multi-FASTA genome
        |
        v
FASTA PARSING
- Read FASTA file
- Extract sequence lengths
- Sort scaffolds in descending order
        |
        v
METRIC CALCULATION
- Calculate Total Assembly Size
- Calculate N50, L50, N90, L90
- Calculate GC content
        |
        v
REPORTING AND OUTPUT
- Generate TXT report
- Generate cumulative bar chart
- Save high-resolution figure (PNG/PDF/SVG)
System Requirements
Minimum Hardware
Component	Recommendation
CPU	2 or more cores recommended
RAM	4 GB minimum; 8 GB or more recommended for large genomes
Storage	Enough space for the input FASTA and output files

Sequences are loaded into memory during processing.



Software Dependencies
Requirement	Details
Operating System	Linux, macOS, or Windows
Python	3.9 or higher
Tkinter	Bundled with most Python distributions
Matplotlib	Required for plotting
NumPy	Required for cumulative calculations


Required Python Modules
Module	Purpose
argparse	Command-line parsing (optional)
os, sys	Filesystem and process utilities
threading	Background execution
tkinter	Graphical user interface
matplotlib	Plotting the cumulative bar chart
numpy	Array and cumulative sum calculations


Installation Guide

Step 1: Verify Python version
bash
python3 --version
Python 3.9 or higher is required.


Step 2: Ensure Tkinter is available
On most systems, Tkinter is bundled with Python. Verify with:

bash
python3 -c "import tkinter; print('Tkinter OK')"
If the import fails on a minimal Linux distribution, install the system package:

bash
sudo apt-get install python3-tk
For Fedora or RHEL:

bash
sudo dnf install python3-tkinter
For Arch:

bash
sudo pacman -S tk


Step 3: Install Matplotlib and NumPy
bash
pip install matplotlib numpy


Step 4: Clone the repository
bash
git clone https://github.com/<your-username>/ assemblyCount.git
cd  assemblyCount


Step 5: Launch the application
bash
python3  assemblyCount.py
Understanding Assembly Metrics
 assemblyCount computes standard genome assembly metrics. The definitions are as follows:



Metric	Definition	Interpretation
- Total Size	Sum of all scaffold lengths.	The total number of base pairs in the assembly.
- N50	The length of the scaffold at which the cumulative sum of scaffolds (from largest to smallest) reaches 50% of the total assembly size.	A measure of contiguity. Higher is better.
- L50	The number of scaffolds required to reach 50% of the total assembly size.	A measure of contiguity. Lower is better.
- N90	The length of the scaffold at which the cumulative sum reaches 90% of the total assembly size.	A measure of completeness. Higher is better.
- L90	The number of scaffolds required to reach 90% of the total assembly size.	A measure of completeness. Lower is better.
- GC Content	The percentage of Guanine and Cytosine bases in the assembly.	Indicates genomic composition and potential contamination.


Input Files
Mandatory Inputs
File	Format	Description
Input FASTA	.fasta, .fa, .fas, .fna	Multi-FASTA genome assembly to analyze


Optional Output Paths
Field	Description
Output Figure	Destination for the generated chart. Defaults to assembly_chart.png
Report file	Destination for the textual report. Defaults to assembly_report.txt



File Preparation
The FASTA must be in standard format; each record starts with the > character.

Only the first whitespace-delimited token after > is used as the scaffold identifier.

Sequences are internally uppercased; case in the input file is not relevant.



Output Files
1. Text Report
The report is a plain text file with the following sections.


Section	Contents
Header	Assembly statistics title
Statistics	Total scaffolds, Total size, N50, L50, N90, L90, GC content
Top Scaffolds	A list of the top 5 largest scaffolds

Example excerpt:

text
--- ASSEMBLY STATISTICS ---
Total Scaffolds: 1542
Total Assembly Size: 450,231,990 bp
N50: 12,450,000 bp (L50: 14)
N90: 2,100,000 bp (L90: 85)
C/G Content: 41.25%

Top 5 Scaffolds (Size in bp):
  1: 45,000,000
  2: 38,200,000
  3: 35,100,000
  4: 32,800,000
  5: 30,500,000



2. High-Resolution Figure
Property	Description
Format	PNG (300 DPI), PDF (Vector), or SVG (Vector)
Chart Type	Cumulative bar chart
Colors	Green (sum of previous scaffolds), Red (current scaffold)
Axes	X-axis: Scaffold Rank (1 = Largest), Y-axis: Cumulative Size (bp)

This figure is designed for publication. The legend, axis labels, and fonts are optimized for readability at high resolutions.


Using the Graphical Interface
Launching the GUI
bash
python3  assemblyCount.py

Main Window Layout
1. Header
Title and one-line description of the tool.


2. Files Panel (Top)
Field	Description
Load FASTA	Browse to select the multi-FASTA to analyze
Process & Save Results	Executes the analysis and prompts for output locations

3. Results Panel (Middle)
Displays the calculated statistics in a monospace text area. It lists the total scaffolds, total size, N50, L50, N90, L90, and GC content, along with the top 5 scaffolds.


4. Plot Panel (Bottom)
Displays the generated cumulative bar chart in real-time after processing.


5. Analysis Log
Real-time console output with color-coded messages.

Green: successful operations

Red: errors

Blue: section headers

Yellow: warnings


Using the Command-Line Interface
While the tool is primarily GUI-based, the core logic can be imported and run in a headless environment. A minimal CLI wrapper is provided for pipeline integration.


Basic Usage
bash
python3  assemblyCount.py --input assembly.fasta --output_chart chart.png --output_report report.txt

Common options
Option	Default	Description
--input	Required	Path to the input FASTA file
--output_chart	assembly_chart.png	Path to save the chart
--output_report	assembly_report.txt	Path to save the text report
--max_scaffolds	50	Number of scaffolds to plot in the chart
Understanding the Results
Interpreting the Chart
The cumulative bar chart is designed to visually represent the composition of the assembly.


The X-axis represents the rank of the scaffold (1 being the largest).

The Y-axis represents the cumulative size in base pairs.

The Green portion of a bar represents the sum of all scaffolds ranked before it.

The Red portion of a bar represents the size of the current scaffold.


If the green portion rises very quickly and the red portions become tiny, it indicates a highly contiguous assembly where a few large scaffolds dominate. If the red portions remain large for many bars, the assembly is fragmented.


Interpreting the Metrics
High N50 / Low L50: Indicates a highly contiguous assembly (good).

Low N50 / High L50: Indicates a fragmented assembly (needs more scaffolding or gap closing).

N90 / L90: Useful for assessing the "long tail" of small scaffolds. A low N90 relative to N50 suggests many small fragments remain.

Troubleshooting
Issue	Solution
Tkinter not found	Install python3-tk (Debian/Ubuntu) or equivalent for your distribution
Matplotlib not found	Run pip install matplotlib
NumPy not found	Run pip install numpy
GUI freezes	Ensure the script is using threading (default in this version). For extremely large genomes (>10GB), consider running from CLI.
Chart looks empty	Ensure the FASTA file contains valid sequences and is not empty.
Out of memory	Extremely large genomes may require more RAM. The tool holds all sequence lengths in memory (not the full sequences), which is usually lightweight.
Chart is too cluttered	The script automatically limits the plot to the top 50 scaffolds. If you need more, change the max_scaffolds variable in the script.



Frequently Asked Questions
Q: What is the difference between N50 and L50?

N50 is a length (e.g., 12 Mb), while L50 is a count (e.g., 14 scaffolds). N50 tells you how long the 50th percentile scaffold is; L50 tells you how many scaffolds it took to reach that 50th percentile. Both are essential for evaluating assembly contiguity.


Q: Does  assemblyCount modify the original FASTA?

No. It only reads the FASTA file to calculate statistics. It does not alter, trim, or edit the original sequences in any way.


Q: Can I generate a chart with more than 50 scaffolds?

Yes, the script limits the chart to the top 50 scaffolds by default to avoid visual clutter. You can modify the max_scaffolds variable in the generate_plot function to increase this number.


Q: What image formats can I save the chart in?

The GUI allows you to save the chart as PNG (300 DPI), PDF (Vector), or SVG (Vector). Vector formats (PDF/SVG) are recommended for publication as they can be scaled infinitely without losing quality.


Q: Does  assemblyCount require external software?

No. It only requires Python 3.9+, Matplotlib, and NumPy. All other modules are part of the standard library.


Q: Can I run  assemblyCount on a server without a display?

Yes. While the GUI requires a display, the core logic can be imported as a library or run via the CLI wrapper to generate the report and chart headlessly.


Q: How long does the analysis take?

For a typical genome with a few thousand scaffolds, the analysis completes in seconds. The runtime scales linearly with the number of scaffolds and the size of the FASTA file.


Citation
If you use  assemblyCount in your research, please cite:

text
 assemblyCount v1.0: A Python-based tool for computing genome assembly statistics. https://github.com/valdirstefenon/assemblyCount


License
 assemblyCount is distributed under the MIT License.

Version History
Version	Date	Changes
v1.0	Sep 2026	Initial public release with GUI, CLI, N50/L50/N90/L90 calculation, GC content, and publication-quality chart generation


assemblyCount v1.0: Making genome assembly assessment accessible to everyone.
